#include "regevEncryption.h"
#include "thresholdEncryption.h"
#include "global.h"
#include "util.h"
#include "seal/seal.h"
#include "seal/util/iterator.h"
#include <numeric>
#include <stdio.h>
#include <NTL/BasicThreadPool.h>
#include <NTL/ZZ.h>
#include <thread>

using namespace seal;
using namespace std;

/*
Each user generates a BFV ciphertexts, one of them encrypting a valid permutation matrix (or part of it, if
ring_dim is not big enough), and all others generate rotations.
When multiplying together, extract each chunk (which are all monomials), multiply, and add them back together.
This will also guarantee that each chunk will eventually also be just monomials.
Each leader is selected by oblivious expanding the monomial in each chunk and multiply with the winner token 
(which will be a ciphertext encryting zero generated based on the public key, this also guarantee the unlink-
ability between public keys and final winner tokens). 
*/

int main() {

    int numcores = 1;
    NTL::SetNumThreads(numcores);
  
    int group_size = 2;

    int committee_size = 2;
    
    // int ring_dim = group_size; // for 200 people, can encrypt ring_dim / 200 Z_p element
    int ring_dim = 32768;
    int n = 4;
    int p = 65537;

    
    int scalar = 1;

    // ideally, we just need the s to be at least 235, set sbar = 400 (400 coefficients, degree of 399), to guarantee sufficient unique s
    // for k = 32, kbar could be 40, set 52 to increase m size and thus reduce repetition time 
    int k = 32, kbar = 40, m = 170, sbar = 300, repetition = 52;
    int coeffToSlot_batch_size = 2*sbar;
    int evaluatePoly_batch_size = sbar;
    int evaluatePoly_party_size = 32768;
    int evaluatePoly_degree = sbar;

    cout << k << kbar << m << sbar << coeffToSlot_batch_size << evaluatePoly_batch_size << evaluatePoly_party_size << evaluatePoly_degree << repetition << scalar << endl;

    EncryptionParameters bfv_params(scheme_type::bfv);
    bfv_params.set_poly_modulus_degree(ring_dim);
    auto coeff_modulus = CoeffModulus::Create(ring_dim, { 
                                                          60, 60, 60, 60, 60,
                                                          60, 60, 60, 60, 60,
                                                          60, 30, 60
                                                        });
    bfv_params.set_coeff_modulus(coeff_modulus);
    bfv_params.set_plain_modulus(p);

    prng_seed_type seed;
    for (auto &i : seed) {
        i = random_uint64();
    }
    auto rng = make_shared<Blake2xbPRNGFactory>(Blake2xbPRNGFactory(seed));
    bfv_params.set_random_generator(rng);

    SEALContext seal_context(bfv_params, true, sec_level_type::none);
    primitive_root = seal_context.first_context_data()->plain_ntt_tables()->get_root();
    // cout << "primitive root: " << primitive_root << endl;

    KeyGenerator new_key_keygen(seal_context, n);
    SecretKey new_key = new_key_keygen.secret_key();
    inverse_ntt_negacyclic_harvey(new_key.data().data(), seal_context.key_context_data()->small_ntt_tables()[0]);

    KeyGenerator keygen(seal_context);
    SecretKey bfv_secret_key = keygen.secret_key();

    MemoryPoolHandle my_pool = MemoryPoolHandle::New();

    PublicKey bfv_public_key;
    keygen.create_public_key(bfv_public_key);

    RelinKeys relin_keys, relin_keys_raise;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(seal_context, bfv_public_key);
    Evaluator evaluator(seal_context);
    BatchEncoder batch_encoder(seal_context);
    Decryptor decryptor(seal_context, bfv_secret_key);

    GaloisKeys gal_keys, gal_keys_expand, gal_keys_coeff, gal_keys_coeff_second;

    vector<int> stepsfirst = {1};
    keygen.create_galois_keys(stepsfirst, gal_keys);

    // gal keys for slotToCoeff
    vector<Modulus> coeff_modulus_coeff = coeff_modulus;
	coeff_modulus_coeff.erase(coeff_modulus_coeff.begin() + 10, coeff_modulus_coeff.end()-1);
	EncryptionParameters parms_coeff = bfv_params;
	parms_coeff.set_coeff_modulus(coeff_modulus_coeff);
	SEALContext context_coeff = SEALContext(parms_coeff, true, sec_level_type::none);

	SecretKey sk_coeff;
	sk_coeff.data().resize(coeff_modulus_coeff.size() * ring_dim);
	sk_coeff.parms_id() = context_coeff.key_parms_id();
	util::set_poly(bfv_secret_key.data().data(), ring_dim, coeff_modulus_coeff.size() - 1, sk_coeff.data().data());
	util::set_poly(
			bfv_secret_key.data().data() + ring_dim * (coeff_modulus.size() - 1), ring_dim, 1,
			sk_coeff.data().data() + ring_dim * (coeff_modulus_coeff.size() - 1));
	KeyGenerator keygen_coeff(context_coeff, sk_coeff);

	vector<int> slotToCoeff_steps_coeff = {0, 1};
	slotToCoeff_steps_coeff.push_back(sqrt(ring_dim/2));
    int offset = coeffToSlot_batch_size;
    while (offset < ring_dim/2) {
        slotToCoeff_steps_coeff.push_back(ring_dim/2 - offset);
        offset += coeffToSlot_batch_size;
    }
    int sq_batch = sqrt(evaluatePoly_batch_size); 
    for (int i = 0; i <= sq_batch; i++) {
        slotToCoeff_steps_coeff.push_back(i * sq_batch);
    }
    slotToCoeff_steps_coeff.push_back(-sq_batch);
    slotToCoeff_steps_coeff.push_back(-evaluatePoly_batch_size);
	keygen.create_galois_keys(slotToCoeff_steps_coeff, gal_keys_coeff);
    keygen.create_relin_keys(relin_keys_raise);


    // gal keys for second slotToCoeff before oblivious expansion
    vector<Modulus> coeff_modulus_coeff_second = coeff_modulus;
	coeff_modulus_coeff_second.erase(coeff_modulus_coeff_second.begin() + 2, coeff_modulus_coeff_second.end()-1);
	EncryptionParameters parms_coeff_second = bfv_params;
	parms_coeff_second.set_coeff_modulus(coeff_modulus_coeff_second);
	SEALContext context_coeff_second = SEALContext(parms_coeff_second, true, sec_level_type::none);

	SecretKey sk_coeff_second;
	sk_coeff_second.data().resize(coeff_modulus_coeff_second.size() * ring_dim);
	sk_coeff_second.parms_id() = context_coeff_second.key_parms_id();
	util::set_poly(bfv_secret_key.data().data(), ring_dim, coeff_modulus_coeff_second.size() - 1, sk_coeff_second.data().data());
	util::set_poly(
			bfv_secret_key.data().data() + ring_dim * (coeff_modulus.size() - 1), ring_dim, 1,
			sk_coeff_second.data().data() + ring_dim * (coeff_modulus_coeff_second.size() - 1));
	KeyGenerator keygen_coeff_second(context_coeff_second, sk_coeff_second);
	keygen_coeff_second.create_galois_keys(slotToCoeff_steps_coeff, gal_keys_coeff_second);


    // gal keys for oblivious expansion
    vector<Modulus> coeff_modulus_expand = coeff_modulus;
    coeff_modulus_expand.erase(coeff_modulus_expand.begin() + 2, coeff_modulus_expand.end()-1);
    EncryptionParameters parms_expand = bfv_params;
    parms_expand.set_coeff_modulus(coeff_modulus_expand);
    SEALContext context_expand = SEALContext(parms_expand, true, sec_level_type::none);

    SecretKey sk_expand;
    sk_expand.data().resize(coeff_modulus_expand.size() * ring_dim);
    sk_expand.parms_id() = context_expand.key_parms_id();
    util::set_poly(bfv_secret_key.data().data(), ring_dim, coeff_modulus_expand.size() - 1, sk_expand.data().data());
    util::set_poly(
            bfv_secret_key.data().data() + ring_dim * (coeff_modulus.size() - 1), ring_dim, 1,
            sk_expand.data().data() + ring_dim * (coeff_modulus_expand.size() - 1));
    KeyGenerator keygen_expand(context_expand, sk_expand);
    vector<uint32_t> galois_elts;
    for (int i = 0; i < ceil(log2(ring_dim)); i++) {
        galois_elts.push_back((ring_dim + exponentiate_uint(2, i)) / exponentiate_uint(2, i));
    }
    keygen_expand.create_galois_keys(galois_elts, gal_keys_expand);


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    Ciphertext random_root_degree_one;
    Plaintext pp;
    pp.resize(ring_dim);
    pp.parms_id() = parms_id_zero;
    vector<uint64_t> perm_v(ring_dim);
    for (int c = 0; c < committee_size; c++) {
        for (int g = 0; g < group_size; g++) {
            for (int i = 0; i < ring_dim; i++) {
                pp.data()[i] = 0;
            }

            pp.data()[1] = 1;
            pp.data()[0] = p - 2*(c+1)*(g+1); // we have {(X-2), (X-4)} for party 1, {(X-4), (X-8)} for party 2


            // batch_encoder.encode(perm_v, pp);
            encryptor.encrypt(pp, random_root_degree_one);
            if (g == 0) cout << "-- [Noise] -- Initial noise: " << decryptor.invariant_noise_budget(random_root_degree_one) << endl;
            saveCiphertext(random_root_degree_one, c * group_size + g);
        }
    }

    // prepare ring_dim different tokens, each is of size ring_dim, in plaintext form
    // for simplicity, use just one token
    Plaintext tokens;
    for (int j = 0; j < (int) ring_dim; j++) {
        perm_v[j] = 2;
    }
    cout << "After param generation." << endl;
    batch_encoder.encode(perm_v, tokens);

    chrono::high_resolution_clock::time_point time_start, time_end;
    uint64_t total_time = 0;

    // into [committee_size] different degree-one poly
    time_start = chrono::high_resolution_clock::now();
   
   
   
   
   
    // vector<Ciphertext> random_root_sum(sbar-1);
    // cout << "Roots: ";

    // // for (int c = 0; c < repetition; c++) {
    //     for (int i = 0; i < sbar-1; i++) {
    //         for (int j = 0; j < ring_dim; j++) {
    //             pp.data()[j] = 0;
    //         }
    //         pp.data()[1] = 1;
    //         pp.data()[0] = (i * 400) % 65537;
    //         cout << (65537 - pp.data()[0]) % 65537 << " ";
    //         encryptor.encrypt(pp, random_root_sum[i]);
    //     }
    // // }
    // cout << endl;

    // for (int i = 0; i < committee_size; i++) {
    //     random_root_sum[i] = EvalAddMany_inpace_modImprove_extract_load(i * group_size, seal_context, bfv_secret_key, committee_size);
    // }

    // time_end = chrono::high_resolution_clock::now();
    // cout << "** [TIME] ** EvalAddMany_inpace_modImprove_extract_load: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    // total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();


    //////////////////////////
    //////////////////////////
    //////////////////////////

    vector<uint64_t> tmp{0,53481,9744,57702,12494,6876,13411,3381,16690,63645,39920,23284,21071,57889,37032,5386,14793,37037,43951,11497,44719,42181,27603,6722,60041,41386,58432,42737,14349,16664,56208,44027,42168,18251,9764,63961,26564,15315,50843,45647,54637,29231,42973,12958,49535,43683,65263,885,16997,51281,18499,11133,62915,59088,14232,50763,31842,1024,41800,30952,46201,34750,38178,47176,28835,8720,38718,25309,11310,48825,6629,48833,60035,22863,18336,63803,58699,11464,40659,10888,13413,43125,26044,36468,16984,1740,58171,68,28165,46747,40733,10691,45820,46940,50149,18966,50993,14079,36990,24446,50627,11364,39108,48254,20725,16139,44237,29999,25034,61587,37368,30081,37616,21402,8708,59506,36646,30809,59899,61415,28982,20110,55850,60976,37588,4827,35505,53054,2535,45206,35662,29517,16882,7999,26942,34477,63825,19820,15049,48777,42449,62581,21323,41076,33705,45741,35605,3258,43086,15147,46543,47595,3630,23146,41212,32288,51857,24980,18278,14733,53294,46232,29573,3611,29852,35623,48875,14318,33427,22890,64916,32769,16269,25390,42035,21518,40651,20800,30538,7997,36583,60253,63995,55754,17126,28108,1531,64379,30357,8715,32231,5479,3982,6198,12042,51541,15777,14135,20511,24971,27675,39288,28896,37131,29937,39567,22912,53488,14003,48398,46153,34311,17202,1347,29301,18857,48536,16411,28014,31815,2636,51143,29070,45161,24599,50100,36122,30403,39662,11112,50976,40986,24018,26701,37994,45058,23246,33621,10341,45825,27996,56662,20959,14071,52987,1038,61574,32254,50901,57483,24824,56455,57252,8036,19564,26373,6584,22274,50392,26402,3286,44503,41534,10306,18405,35175,1987,5169,54374,27025,54228,16059,55638,53709,64557,31144,26258,61614,38522,43658,56775,15860,14495,20001,58060,61697,54055,1701,9693,31755,5139,30306,6881,28543,35258,13643,10831,50589,45442,33444,24198,11185,52602,26574,48960,49279,27671,31821,7525,59322,8331,65491,46888,16466,5783,26968,27500,16060,56366,38910,12731,22345,33194,8317,20913,46616,17843,4053,52027,29151,21178,52685,23430,42051,32755,61419,39115,60340,63508,18229,51439,54229,53588,28016,15294,39158,46981,28791,63620,16525,24938,19385,41615,23334,18274,22209,6898,58017,57534,37310,58294,49219,22106,45492,24435,40315,52358,1850,13879,29983,57676,58029,63140,10683,34125,19300,28435,59951,17253,56215,14411,37098,39457,18950,65083,11951,39444,37739,61958,57502,64858,50677,25197,34573,11418,41859,41883,16093,40492,1};
    vector<uint64_t> vvv(ring_dim);
    for (int i = 0; i < ring_dim/2; i++) {
        if (i >= 16000) {
            vvv[i] = 0;
        } else if ((i / 400) % 2 == 0) {
            vvv[i] = tmp[i % 400];
        } else {
            vvv[i] = 0;
        }
        vvv[i + ring_dim/2] = vvv[i];
    }
    batch_encoder.encode(vvv,pp);
    vector<Ciphertext> random_poly(1);
    encryptor.encrypt(pp, random_poly[0]);

    // print_ct_to_vec(random_poly[0], seal_context, bfv_secret_key, ring_dim);

    cout << decryptor.invariant_noise_budget(random_poly[0]) << endl;

    vector<vector<Ciphertext>> results = evaluatePolynomial_batch(bfv_secret_key, seal_context, random_poly, gal_keys_coeff, ring_dim,
                                                                  evaluatePoly_batch_size, evaluatePoly_party_size, evaluatePoly_degree,
                                                                  numcores);

    cout << decryptor.invariant_noise_budget(results[0][0]) << endl;
    //////////////////////////
    //////////////////////////
    //////////////////////////

    // // into one single poly with [committee_size] different roots
    // Ciphertext final_random_root;
    // time_start = chrono::high_resolution_clock::now();
    // NTL::SetNumThreads(8);
    // NTL_EXEC_RANGE(1, first, last);
    // for (int c = first; c < last; c++) {
    //     final_random_root = EvalMultMany_inpace_modImprove(random_root_sum, relin_keys, seal_context, bfv_secret_key);
    // }
    // NTL_EXEC_RANGE_END;

    // print_ct_to_pl(final_random_root, seal_context, bfv_secret_key, ring_dim);



    // time_end = chrono::high_resolution_clock::now();
    // cout << "** [TIME] ** EvalMultMany_inpace_modImprove: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    // total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

    // cout << "-- [Noise] -- After: " << decryptor.invariant_noise_budget(final_random_root) << endl;

    // evaluator.mod_switch_to_next_inplace(final_random_root);
    // cout << "-- [Noise] -- After +++ check mod: " << decryptor.invariant_noise_budget(final_random_root) << endl;

    // Ciphertext random_poly = final_random_root;
    // // evaluate the poly to extract all roots on indices
    // time_start = chrono::high_resolution_clock::now();
    // vector<Ciphertext> randoms(repetition);
    // for (int i = 0; i < repetition; i++) {
    //     randoms[i] = random_poly;
    // }



    // vector<Ciphertext> random_poly_on_slots = coeffToSlot_WOPreprocess_batch(bfv_secret_key, seal_context, randoms, gal_keys_coeff,
    //                                                                           ring_dim, coeffToSlot_batch_size, p, scalar, numcores);
    // time_end = chrono::high_resolution_clock::now();
    // cout << "** [TIME] ** coeffToSlot_WOPrepreocess: " << 
    //       chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    // total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

    // // print_ct_to_vec(random_poly_on_slots[0], seal_context, bfv_secret_key, ring_dim);







    
    // cout << "-- [Noise] -- After slotToCoeff: " << decryptor.invariant_noise_budget(random_poly_on_slots[0]) << endl;
    // for (int i = 0; i < (int) random_poly_on_slots.size(); i++) {
    //     evaluator.mod_switch_to_next_inplace(random_poly_on_slots[i]);
    // }
    // cout << "-- [Noise] -- After +++ check mod: " << decryptor.invariant_noise_budget(random_poly_on_slots[0]) << endl;


    // vector<vector<Ciphertext>> binary_poly_roots = evaluatePolynomial_batch(bfv_secret_key, seal_context, random_poly_on_slots, gal_keys_coeff, ring_dim,
    //                                                                 evaluatePoly_batch_size, evaluatePoly_party_size, evaluatePoly_degree);



























    // time_start = chrono::high_resolution_clock::now();
    // Ciphertext binary_vector = evaluatePolynomial(seal_context, random_poly, gal_keys_coeff, ring_dim, 8);
    // time_end = chrono::high_resolution_clock::now();
    // cout << "** [TIME] ** evaluatePolynomial: " <<
    //       chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    // total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

    // cout << "-- [Noise] -- After evaluate poly: " << decryptor.invariant_noise_budget(binary_vector) << endl;

    // evaluator.mod_switch_to_next_inplace(binary_vector);

    // cout << "-- [Noise] -- After evaluate poly +++ check mod: " << decryptor.invariant_noise_budget(binary_vector) << endl;


    map<int, bool> raise_mod = {{4, false}, {16, false}, {64, false}, {256, false}};
    time_start = chrono::high_resolution_clock::now();
    Ciphertext binary_vector = raisePowerToPrime(seal_context, relin_keys_raise, results[0][0], raise_mod, raise_mod, 256, 256, p); // all one except v-th slot is zero
    time_end = chrono::high_resolution_clock::now();
    cout << "** [TIME] ** raisePowerToPrime: " <<
          chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    // total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();


    // // flip binary bits
    // vector<uint64_t> intInd(ring_dim, 1); 
	// Plaintext plainInd;
	// batch_encoder.encode(intInd, plainInd);
	// evaluator.negate_inplace(binary_vector);
	// evaluator.add_plain_inplace(binary_vector, plainInd);

    // time_start = chrono::high_resolution_clock::now();
    // binary_vector = slotToCoeff_WOPrepreocess(seal_context, binary_vector, gal_keys_coeff, ring_dim, p, 1, 1);
    // time_end = chrono::high_resolution_clock::now();
    // cout << "** [TIME] ** slotToCoeff_WOPrepreocess: " <<
    //       chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    // total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

    // print_ct_to_pl(binary_vector, seal_context, bfv_secret_key, ring_dim);

    // cout << "now the noise? " << decryptor.invariant_noise_budget(binary_vector) << endl;

    cout << "\n\nTotal time: " << total_time << endl;


    // vector<vector<uint64_t>> U = generateMatrixU_transpose(4, p);
    // vector<vector<uint64_t>> U_inverse = generateInverse_vander(U, p);

    // cout << U << endl << U_inverse << endl;



    return 0;
}
