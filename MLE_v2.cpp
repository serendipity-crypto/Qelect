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

    int numcores = 8;
    NTL::SetNumThreads(numcores);
  
    int group_size = 1024;

    // int committee_size = 2;
    
    // int ring_dim = group_size; // for 200 people, can encrypt ring_dim / 200 Z_p element
    int ring_dim = 32768;
    int n = 32768;
    int p = 65537;

    
    int scalar = 1;

    EncryptionParameters bfv_params(scheme_type::bfv);
    bfv_params.set_poly_modulus_degree(ring_dim);
    auto coeff_modulus = CoeffModulus::Create(ring_dim, { 
                                                          60, 60, 60, 60,
                                                          60, 60, 60, 60,
                                                          60, 60, 60, 60
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

    vector<int> stepsfirst = {0, 1};
    for (int i = 0; i < log2(ring_dim/2); i++) {
        stepsfirst.push_back(-(1<<i));
        stepsfirst.push_back((1<<i));
    }
    keygen.create_galois_keys(stepsfirst, gal_keys);

    // gal keys for slotToCoeff
    vector<Modulus> coeff_modulus_coeff = coeff_modulus;
	// coeff_modulus_coeff.erase(coeff_modulus_coeff.begin() + 13, coeff_modulus_coeff.end()-1);
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
	keygen_coeff.create_galois_keys(slotToCoeff_steps_coeff, gal_keys_coeff);
    keygen_coeff.create_relin_keys(relin_keys_raise);


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

    // uint64_t total_time = 0;


    // prepare ring_dim different tokens, each is of size ring_dim, in plaintext form
    // for simplicity, use just one token
    Plaintext tokens;
    vector<uint64_t> vvv(ring_dim);
    for (int j = 0; j < (int) ring_dim; j++) {
        vvv[j] = random_uint64() % p;
    }
    cout << "After param generation." << endl;
    batch_encoder.encode(vvv, tokens);

    chrono::high_resolution_clock::time_point time_start, time_end;    

	Ciphertext tmp_random;
    Plaintext pp;






    // for (int g = 0; g < group_size; g++) {
    //     for (int i = 0; i < ring_dim; i++) {
    //         vvv[i] = random_uint64() % p;
    //     }
    //     batch_encoder.encode(vvv, pp);
    //     encryptor.encrypt(pp, tmp_random);
    //     saveCiphertext(tmp_random, g);
    // }

    // cout << "After preparing ciphertexts...\n";
    // cout << "Initial noise: " << decryptor.invariant_noise_budget(tmp_random) << endl;


    // time_start = chrono::high_resolution_clock::now();
    // NTL::SetNumThreads(numcores);
    // vector<Ciphertext> subsum_random_core(numcores);
    // uint64_t batch_add_size = group_size / numcores;
    // NTL_EXEC_RANGE(numcores, first, last);
    // for (int c = first; c < last; c++) {
    //     subsum_random_core[c] = EvalAddMany_inpace_modImprove_extract_load(c * batch_add_size, seal_context, bfv_secret_key, batch_add_size);
    // }
    // NTL_EXEC_RANGE_END;

    // time_end = chrono::high_resolution_clock::now();
    // total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    // cout << "       after sum up: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;

    // Ciphertext random = EvalAddMany_inpace_modImprove_extract_multi_core(subsum_random_core, seal_context, bfv_secret_key);
    // cout << "Noise after sum up: " << decryptor.invariant_noise_budget(random) << endl;
    for (int i = 0; i < ring_dim; i++) {
        vvv[i] = random_uint64() % group_size;
    }
    batch_encoder.encode(vvv, pp);
    encryptor.encrypt(pp, tmp_random);
    Ciphertext random = tmp_random;




    
    time_start = chrono::high_resolution_clock::now();
	Ciphertext extracted;
   	for (int i = 0; i < ring_dim; i++) {
		vvv[i] = 0;
    }
	vvv[0] = 1; // extractor
	batch_encoder.encode(vvv, pp);
	evaluator.multiply_plain(random, pp, extracted);

    cout << "Tested random indices... : ";
    print_ct_to_vec(random, seal_context, bfv_secret_key, 10);
    
    cout << "Initial noise budget: " << decryptor.invariant_noise_budget(extracted) << endl;

	// expand random values to all slots
    expandFirstToAll(bfv_secret_key, seal_context, extracted, gal_keys, ring_dim);
    Ciphertext extracted_all_slots = extracted;
    // print_ct_to_vec(extracted, seal_context, bfv_secret_key, 1000);

	// Ciphertext extracted_all_slots = slotToCoeff_WOPrepreocess(seal_context, extracted, gal_keys_coeff, ring_dim, p, scalar, numcores);
	time_end = chrono::high_resolution_clock::now();
    total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
	// cout << "	extract and expand to all dimension: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
	// cout << "after fill: " << decryptor.invariant_noise_budget(extracted_all_slots) << endl;
	
    time_start = chrono::high_resolution_clock::now();
	evaluator.negate_inplace(extracted_all_slots); // [-r, -r, ..., -r]
    map<int, bool> raise_mod_1 = {{2, false}, {8, false}, {32, false}, {128, false}};
	map<int, bool> raise_mod_2 = {{4, false}, {16, false}, {64, false}, {256, false}};

    Plaintext all_ones;
    all_ones.resize(ring_dim);
    all_ones.parms_id() = parms_id_zero;
    for (int i = 0; i < ring_dim; i++) {
        all_ones.data()[i] = 0;
    }
    all_ones.data()[0] = 1;

	// vector<Ciphertext> evaluated(2);
    // NTL_EXEC_RANGE(2, first, last); // ring_dim = 32768, evaluate points: 65536, two batches
    // for (int c = first; c < last; c++) {
    //     for (int i = 0; i < ring_dim; i++) {
    //         vvv[i] = i + ring_dim*c;
    //     }
    //     batch_encoder.encode(vvv, pp);
    //     evaluator.add_plain(extracted_all_slots, pp, evaluated[c]); // [-r, -r, ..., -r] + [0, 1, ..., n-1]
    //     evaluated[c] = raisePowerToPrime(seal_context, relin_keys, evaluated[c], raise_mod_2, raise_mod_2, 256, 256, p);
    //     evaluator.negate_inplace(evaluated[c]);
    //     evaluator.add_plain_inplace(evaluated[c], all_ones);
    // }
    // NTL_EXEC_RANGE_END;

    Ciphertext evaluated;
    for (int i = 0; i < ring_dim; i++) {
        vvv[i] = i;
    }
    batch_encoder.encode(vvv, pp);
    evaluator.add_plain(extracted_all_slots, pp, evaluated); // [-r, -r, ..., -r] + [0, 1, ..., n-1]
    evaluated = raisePowerToPrime(seal_context, relin_keys, evaluated, raise_mod_2, raise_mod_2, 256, 256, p);
    evaluator.negate_inplace(evaluated);
    evaluator.add_plain_inplace(evaluated, all_ones);

	time_end = chrono::high_resolution_clock::now();
    total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
	// cout << "	raisePower: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
	// cout << decryptor.invariant_noise_budget(evaluated) << endl;

    time_start = chrono::high_resolution_clock::now();
    sumUpEvaluation_ToPartyRange(seal_context, evaluated, group_size, gal_keys);
    time_end = chrono::high_resolution_clock::now();
    total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
	// cout << "	sumUpEvaluation_ToPartyRange: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;

    // evaluator.mod_switch_to_next_inplace(evaluated);
	
    // NOTICE!!!!!!!!!!!!!! WE COULD DO BATCHING HERE!!!!!!!!!!!!!!!!!!!!! YEAH!!!!!!!!!!!!!!!!!!!!!

    // time_start = chrono::high_resolution_clock::now();
    evaluated = slotToCoeff_WOPrepreocess(seal_context, evaluated, gal_keys_coeff, ring_dim, p, scalar, numcores, group_size);
	// time_end = chrono::high_resolution_clock::now();
    // total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
	// cout << "	second: slotToCoeff_WOPrepreocess: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
	// cout << decryptor.invariant_noise_budget(evaluated) << endl;

	// evaluator.mod_switch_to_next_inplace(evaluated);
    evaluator.mod_switch_to_next_inplace(evaluated);
	// cout << evaluated.coeff_modulus_size() << endl;
	// cout << decryptor.invariant_noise_budget(evaluated) << endl;
    // print_ct_to_pl(evaluated, seal_context, bfv_secret_key, 128);

    time_start = chrono::high_resolution_clock::now();
    vector<Ciphertext> expanded_for_batch(ring_dim/group_size); // how many batch one ring can hold, and each time we only use one of them

    if (ring_dim/group_size == 1) {
        expanded_for_batch[0] = evaluated;
    } else {
        expanded_for_batch = subExpand(sk_expand, context_expand, parms_expand, evaluated, ring_dim, gal_keys_expand,
                                       ring_dim/group_size);
    }

    time_end = chrono::high_resolution_clock::now();
    total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() / (ring_dim/group_size);
    // cout << "       First level subExpand time: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;

    // ready for multi-core acceleration
    time_start = chrono::high_resolution_clock::now();
    vector<Ciphertext> expanded_subtree_roots = subExpand(sk_expand, context_expand, parms_expand, expanded_for_batch[0], ring_dim,
                                                                 gal_keys_expand, numcores);
    time_end = chrono::high_resolution_clock::now();
    total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    // cout << "       second expand core time: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;

    vector<Ciphertext> cached(numcores);
    for (int i = 0; i < (int) expanded_subtree_roots.size(); i++) {
        // evaluator.mod_switch_to_next_inplace(expanded_subtree_roots[i]);
        cached[i] = expanded_subtree_roots[i];
    }

    // print_ct_to_vec(expanded_subtree_roots[0], seal_context, bfv_secret_key, 100);

    time_start = chrono::high_resolution_clock::now();
    vector<vector<Ciphertext>> expanded_leaf(numcores, vector<Ciphertext>(group_size/numcores));  
    vector<Ciphertext> token_subsum(numcores);                                             
    NTL_EXEC_RANGE(numcores, first, last);
    for (int i = first; i < last; i++) {
        // cout << i << endl; 
        // expand each 1 out of the 8 subroots to leaf level
        expanded_leaf[i] = expand(context_expand, parms_expand, expanded_subtree_roots[i], ring_dim,
                                  gal_keys_expand, group_size/numcores);

        for (int j = 0; j < (int) expanded_leaf[i].size(); j++) {
            // evaluator.multiply_plain(tmp_random, tokens, expanded_leaf[i][j]);
            evaluator.multiply_plain_inplace(expanded_leaf[i][j], tokens);
        }
        token_subsum[i] = EvalAddMany_inpace_modImprove_extract_multi_core(expanded_leaf[i], seal_context, bfv_secret_key);
    }
    NTL_EXEC_RANGE_END;

    Ciphertext result = EvalAddMany_inpace_modImprove_extract_multi_core(token_subsum, seal_context, bfv_secret_key);
    // print_ct_to_vec(result, seal_context, bfv_secret_key, 100);

    vector<Ciphertext> final_par_dec(group_size);
    for (int i = 0; i < group_size; i++) {
        final_par_dec[i] = result;
    }

    result = EvalAddMany_inpace_modImprove_extract_multi_core(token_subsum, seal_context, bfv_secret_key);

    stringstream data_streamdg;
    // cout << result.coeff_modulus_size() << endl;
    auto digsize = result.save(data_streamdg);

	// simulate the file saving
	ofstream final_ct_bytes, part_dec_bytes;

	final_ct_bytes.open("../data/final_ct_bytes.txt");
	final_ct_bytes << data_streamdg.rdbuf();
	final_ct_bytes.close();

	result.save(data_streamdg)
	part_dec_bytes.open("../data/final_ct_bytes.txt");
	part_dec_bytes << data_streamdg.rdbuf();
	part_dec_bytes.close();
	


    time_end = chrono::high_resolution_clock::now();
    total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    // cout << "** [TIME] ** Second + Leaf level expand + multi token time: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    // total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();


    // vector<Ciphertext> expanded_leaf = subExpand(sk_expand, context_expand, parms_expand, evaluated[0], ring_dim, gal_keys_expand, ring_dim/2);
    // cout << "After this?\n";
    // // vector<Ciphertext> expanded_leaf = expand(context_expand, parms_expand, expandeddd[0], ring_dim, gal_keys_expand, ring_dim/4);

	// time_end = chrono::high_resolution_clock::now();
	// cout << "	expand: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
	// cout << decryptor.invariant_noise_budget(expanded_leaf[0]) << endl;

	// Ciphertext result;
	// for (int i = 0; i < (int) expanded_leaf.size(); i++) {
	// 	if (i == 0) {
	// 		evaluator.multiply_plain(expanded_leaf[i], pp, result);
	// 	} else {
	// 		evaluator.multiply_plain_inplace(expanded_leaf[i], pp);
	// 		evaluator.add_inplace(result, expanded_leaf[i]);
	// 	}
	// }
	cout << "Lefted noise budget: " << decryptor.invariant_noise_budget(result) << endl;
    cout << "Total broadcast communication size (for each party): " << digsize * group_size * 2 / 1000000 << " MB." << endl;

    cout << "\n---------------------------------\n";
    cout << "Total number of parties: " << group_size << endl;
    cout << "Preprocessed time      : " << preprocess_time << " us." << endl;
	cout << "Total time             : " << total_time << " us." << endl;

    return 0;
}
