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

    // ideally, we just need the s to be at least 235, set sbar = 400, to guarantee sufficient unique s
    // for k = 32, kbar could be 40, set 52 to increase m size and thus reduce repetition time 
    int k = 32, kbar = 52, m = 210, sbar = 400, repetition = 40;
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
   
   
   
   
   
    vector<Ciphertext> random_root_sum(sbar);
    // cout << "Roots: ";

    // for (int c = 0; c < repetition; c++) {
        for (int i = 0; i < sbar; i++) {
            for (int j = 0; j < ring_dim; j++) {
                pp.data()[j] = 0;
            }
            pp.data()[1] = 1;
            pp.data()[0] = i * 80;
            encryptor.encrypt(pp, random_root_sum[i]);
        }
    // }

    // for (int i = 0; i < committee_size; i++) {
    //     random_root_sum[i] = EvalAddMany_inpace_modImprove_extract_load(i * group_size, seal_context, bfv_secret_key, committee_size);
    // }

    // time_end = chrono::high_resolution_clock::now();
    // cout << "** [TIME] ** EvalAddMany_inpace_modImprove_extract_load: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    // total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();



    vector<uint64_t> tmp{0,14220,39374,46529,18277,40538,48595,35956,24748,29768,19368,29599,40654,9422,39143,63366,5374,39868,42622,10381,4701,8564,30189,17032,62351,30577,17106,43456,49003,18972,42549,1958,46908,36074,46194,21772,27831,48054,50993,58641,16990,26484,25834,62978,14065,12282,5056,46016,7133,46858,39451,54544,46128,13962,9465,10533,36411,1133,65026,26807,33268,22314,49146,25707,64838,24129,57954,62167,16802,43421,25729,17633,64013,62073,426,7582,46885,10962,47666,34424,13077,3577,61511,8565,31953,14352,28872,1113,32163,62390,50476,46576,22794,38260,30806,18847,40624,51017,5374,3777,48521,51281,40398,21886,39288,34548,7088,40840,47704,27118,8220,63239,64925,64287,38199,38141,25519,24386,19320,53291,18399,8075,14310,28144,20804,16243,44064,62707,47937,26454,63627,17999,23954,4049,7478,14425,12274,38244,42219,22953,50802,7600,14865,31288,46072,28780,59118,35701,36414,1909,8508,37949,44007,2075,14103,13705,8137,26243,54843,54466,32490,35896,4656,54019,20721,30498,39064,43275,46063,58076,37476,39000,14805,15811,22502,43312,45793,9969,11680,22781,56469,6716,56224,55162,26495,25034,13148,13129,24151,63177,52001,2226,29507,21154,51498,15088,33876,52468,52354,32214,10123,5797,65435,54017,64588,46244,16702,52744,14370,58434,29813,33750,19756,28071,40144,6398,32425,35753,29882,33726,52259,29677,42005,46085,52826,7703,64499,44817,27471,36675,63764,53826,20882,18240,679,23446,42036,41739,62235,14676,8671,64770,26802,32989,43723,890,4355,60391,13839,63778,52925,23963,3072,1827,15007,50969,4031,19583,26924,27747,7113,33489,26797,15453,40939,21894,5785,47334,57852,8231,14819,12762,3827,6910,30465,3234,28106,14417,41147,25947,60269,19143,48801,56552,58138,356,60882,42681,39361,58853,27275,44896,23880,1552,23325,21864,65503,46955,2964,47906,10166,19783,141,31359,56768,11321,26408,18977,43997,60420,39881,19178,42009,2810,40623,29467,754,7218,4586,65112,47955,34052,57889,22298,41663,15874,63394,54792,51167,364,17533,39516,20181,11932,34698,13661,1867,17012,21285,27716,24242,34656,43952,61316,41986,19195,17673,18582,63438,44478,10622,19464,27882,13611,4225,47395,8295,63596,3052,13725,53749,23728,57135,19159,38643,59434,11513,30737,29849,40350,25974,20443,22194,10603,53674,15098,457,52900,53987,5116,58202,23293,31365,55770,24133,35897,28705,8304,41283,17125,11708,5431,55792,26520,42478,52932,56610,3863,52480,26911};
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

    vector<vector<Ciphertext>> results = evaluatePolynomial_batch(bfv_secret_key, seal_context, random_poly, gal_keys_coeff, ring_dim,
                                                                  evaluatePoly_batch_size, evaluatePoly_party_size, evaluatePoly_degree);


    // into one single poly with [committee_size] different roots
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


    // map<int, bool> raise_mod = {{4, false}, {16, false}, {64, false}, {256, false}};
    // time_start = chrono::high_resolution_clock::now();
    // binary_vector = raisePowerToPrime(seal_context, relin_keys_raise, binary_vector, raise_mod, raise_mod, 256, 256, p); // all one except v-th slot is zero
    // time_end = chrono::high_resolution_clock::now();
    // cout << "** [TIME] ** raisePowerToPrime: " <<
    //       chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
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
