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
   
   
   
   
   
    // vector<Ciphertext> random_root_sum(sbar-1);
    // // cout << "Roots: ";

    // // for (int c = 0; c < repetition; c++) {
    //     for (int i = 0; i < sbar-1; i++) {
    //         for (int j = 0; j < ring_dim; j++) {
    //             pp.data()[j] = 0;
    //         }
    //         pp.data()[1] = 1;
    //         pp.data()[0] = (i * 100) % 65537;
    //         encryptor.encrypt(pp, random_root_sum[i]);
    //     }
    // // }

    // for (int i = 0; i < committee_size; i++) {
    //     random_root_sum[i] = EvalAddMany_inpace_modImprove_extract_load(i * group_size, seal_context, bfv_secret_key, committee_size);
    // }

    // time_end = chrono::high_resolution_clock::now();
    // cout << "** [TIME] ** EvalAddMany_inpace_modImprove_extract_load: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    // total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();


    //////////////////////////
    //////////////////////////
    //////////////////////////

    vector<uint64_t> tmp{0,3715,33783,25887,14141,48723,45800,62156,64314,30272,1063,3163,50506,65059,9258,5386,59172,2759,60310,59604,47430,17844,42252,58815,21984,58731,61498,4007,52449,33810,14052,44027,37598,29868,35063,55303,3681,11331,36442,19890,43600,56600,2282,25139,1798,55979,32700,885,2451,34052,4270,31957,2089,61844,61979,14774,3706,49153,11817,6265,7890,10364,42313,47176,49803,8446,53083,56478,46928,33813,14727,16704,22008,27414,6162,50682,55190,33485,26549,10888,53652,34630,28391,29554,24311,49044,34610,65469,18414,38492,14568,15658,4812,19318,61690,18966,7361,28653,8028,32161,2281,15674,55760,17283,48174,3924,52460,53622,55688,57098,9342,30081,19390,14747,33016,28952,38340,34939,34178,4122,15146,5925,30135,53487,45644,20782,58029,53054,10140,2389,54110,19597,50937,60941,26033,31060,6848,10565,19919,30655,48792,48968,21715,41076,3746,10949,50462,47604,13663,44110,37517,17942,51017,22886,49449,57471,48939,50714,37338,14733,16565,18805,57636,6898,28206,26446,36934,51219,62903,26982,39744,65409,52479,9779,26893,21518,31530,5115,53859,15585,39365,49483,33154,9783,62570,9031,33090,34300,44507,21025,24442,5479,15928,33631,49781,21559,33546,27789,44025,40566,20374,26762,51229,62866,15828,6569,5728,53488,56012,53461,4627,1658,50932,12204,9059,46680,2467,65109,42140,47485,53290,40061,40036,45161,32859,15156,18013,49802,46485,32074,52793,24551,35002,31543,58790,65201,51564,47158,51738,45825,46447,54611,30636,63178,59589,57280,17375,33283,58544,63327,49689,31197,29567,49655,4891,26373,26336,28699,13775,8601,22477,25891,22385,55231,57454,27033,3906,53013,27474,63130,13557,16059,25941,7363,2817,42887,18022,53494,23138,21879,35048,8388,55375,57167,54156,65297,29898,1701,38772,49321,1211,24970,33685,59657,23954,51894,22213,42557,40877,23683,59771,62140,45919,26574,64766,2020,1445,19588,37771,37253,47070,46,9059,64229,23110,43114,20910,17388,46860,38910,50924,29835,27232,31968,49850,29855,44692,61484,54040,57880,20885,13262,59759,55877,24573,61419,25386,47922,1218,13497,47325,17091,52140,37521,4361,28842,7918,35185,62435,13321,39003,19385,35386,45659,55407,49322,51093,470,18385,28227,28972,64477,27030,19634,13694,23000,45858,1850,55516,20969,21192,44062,35878,44389,7853,46237,17334,23839,9937,27100,54498,26895,59017,18950,63721,60142,34010,27245,5176,53751,16554,14860,30286,36665,55692,32164,38543,13294,10123,1};
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
