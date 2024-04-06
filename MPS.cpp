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


int main() {
  
    int ring_dim = 8192; // for 200 people, can encrypt ring_dim / 200 Z_p element
    int n = 512;
    int p = 65537;

	int numcores = 2;

    int group_size = 256;
    int batch_size = ring_dim / group_size;
    // int batch_size = ring_dim / group_size;

    EncryptionParameters bfv_params(scheme_type::bfv);
    bfv_params.set_poly_modulus_degree(ring_dim);
    auto coeff_modulus = CoeffModulus::Create(ring_dim, { 30, 60, 60, 60,  
                                                          60, 60, 60, 60,
                                                          60 });
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
    cout << "primitive root: " << primitive_root << endl;

    KeyGenerator new_key_keygen(seal_context, n);
    SecretKey new_key = new_key_keygen.secret_key();
    inverse_ntt_negacyclic_harvey(new_key.data().data(), seal_context.key_context_data()->small_ntt_tables()[0]);

    KeyGenerator keygen(seal_context);
    SecretKey bfv_secret_key = keygen.secret_key();

    MemoryPoolHandle my_pool = MemoryPoolHandle::New();

    // generate a key switching key based on key_before and secret_key
    KSwitchKeys ksk;
    seal::util::ConstPolyIter secret_key_before(bfv_secret_key.data().data(), ring_dim, coeff_modulus.size());

    new_key_keygen.generate_kswitch_keys(secret_key_before, 1, static_cast<KSwitchKeys &>(ksk), false);
    ksk.parms_id() = seal_context.key_parms_id();

    PublicKey bfv_public_key;
    keygen.create_public_key(bfv_public_key);

    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(seal_context, bfv_public_key);
    Evaluator evaluator(seal_context);
    BatchEncoder batch_encoder(seal_context);
    Decryptor decryptor(seal_context, bfv_secret_key);

    GaloisKeys gal_keys;
    vector<int> rot_steps = {1};
    for (int i = 0; i < n;) {
        rot_steps.push_back(i);
        i += sqrt(n);
    }
    for (int i = 1; i < batch_size; i++) {
      rot_steps.push_back(i);
    }
    for (int i = 0; i < ring_dim/2;) {
        if (find(rot_steps.begin(), rot_steps.end(), i) == rot_steps.end()) {
            rot_steps.push_back(i);
        }
        i += sqrt(ring_dim/2);
    }
    keygen.create_galois_keys(rot_steps, gal_keys);


    GaloisKeys gal_keys_slotToCoeff;
    vector<int> slotToCoeff_steps_coeff = {0, 1};
    slotToCoeff_steps_coeff.push_back(sqrt(ring_dim/2));
    keygen.create_galois_keys(slotToCoeff_steps_coeff, gal_keys_slotToCoeff);


    GaloisKeys gal_keys_expand;
    vector<uint32_t> galois_elts;
    for (int i = 0; i < ceil(log2(ring_dim)); i++) {
      galois_elts.push_back((ring_dim + exponentiate_uint(2, i)) / exponentiate_uint(2, i));
    }
    keygen.create_galois_keys(galois_elts, gal_keys_expand);

    // // [0,0,1,0] * [0,1,0,0]
    // Plaintext plainInd1, plainInd2;
    // plainInd1.resize(ring_dim);
    // plainInd1.parms_id() = parms_id_zero;
    // for (int i = 0; i < (int) ring_dim; i++) {
    //     plainInd1.data()[i] = 0;
    // }
    // plainInd1.data()[2] = 1;

    // plainInd2.resize(ring_dim);
    // plainInd2.parms_id() = parms_id_zero;
    // for (int i = 0; i < (int) ring_dim; i++) {
    //     plainInd2.data()[i] = 0;
    // }
    // plainInd2.data()[3] = 65536;

    // Ciphertext c;
    // encryptor.encrypt(plainInd2, c);
    // evaluator.multiply_plain_inplace(c, plainInd1);

    // decryptor.decrypt(c, plainInd2);
    // for (int i = 0; i < (int) ring_dim; i++) {
    //     cout << plainInd2.data()[i] << " ";
    // }
    // cout << endl;


    // cout << "-----------------------------------------------------\n";



    // Plaintext ppp;
    // ppp.resize(ring_dim);
    // ppp.parms_id() = parms_id_zero;
    // for (int j = 0; j < (int) ring_dim; j++) {
    //   ppp.data()[j] = 1;
    // }
    // ppp.data()[101] = 1;
    // ppp.data()[10] = 1;
    // Ciphertext ccc;
    // encryptor.encrypt(ppp, ccc);
    

    // map<int, bool> raise_mod = {{4, false}, {16, false}, {64, false}, {256, false}};
    // Ciphertext tmp = raisePowerToPrime(seal_context, relin_keys, ccc, raise_mod, raise_mod,
    // 				       256, 256, p);

    // decryptor.decrypt(tmp, ppp);
    // for (int j = 0; j < (int) ring_dim; j++) {
    //   cout << ppp.data()[j] << " ";
    // }
    // cout << endl;


    // cout << "-----------------------------------------------------\n";




    
    // prepare 256 permutation vectors, representing the first row of a permutation matrix
    vector<Ciphertext> perms(group_size);

    Plaintext pp;
    pp.resize(ring_dim);
    pp.parms_id() = parms_id_zero;
    for (int i = 0; i < (int) perms.size(); i++) {
		for (int j = 0; j < (int) ring_dim; j++) {
			pp.data()[j] = 0;
		}
		if (i == 0) {
			pp.data()[2] = 1; // [0, 0, 1, 0,...]
			// for (int b = 0; b < batch_size; b++) {
			//   pp.data()[2 + b * group_size] = 1; // [0,0,1,0]
			// }
		} else {
			pp.data()[0] = 1; // [1,0,0,0], notice that following ciphertext no need to repeat
		}
		encryptor.encrypt(pp, perms[i]);
    }

    cout << "--[NOISE]-- initial: " << decryptor.invariant_noise_budget(perms[0]) << endl;



    // prepare ring_dim different tokens, each is of size ring_dim, in plaintext form
    vector<Plaintext> tokens(ring_dim); 
    for (int i = 0; i < (int) tokens.size(); i++) {
		tokens[i].resize(ring_dim);
		tokens[i].parms_id() = parms_id_zero;
		for (int j = 0; j < (int) ring_dim; j++) {
			tokens[i].data()[j] = random_uint64() % 65537;
			// if (i == 2) {
			//   cout << tokens[i].data()[j] << " ";
			// }
		}
    }
    cout << "After generation" << endl;

    chrono::high_resolution_clock::time_point time_start, time_end;
	uint64_t total_time = 0;
    time_start = chrono::high_resolution_clock::now();

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // multiply all perm vec ciphertexts, log(n) depth of ct-multi

    NTL::SetNumThreads(numcores);
	vector<Ciphertext> perm_share_final(numcores);
	uint64_t batch_perm_size = perms.size() / numcores;
    NTL_EXEC_RANGE(numcores, first, last);
	vector<Ciphertext>::const_iterator b = perms.begin();
    for (int i = first; i < last; i++) {
		vector<Ciphertext> perms_share(b + i*batch_perm_size , b + (i+1)*batch_perm_size);
    	perm_share_final[i] = EvalMultMany_inpace_modImprove(perms_share, relin_keys, seal_context, bfv_secret_key);
    }
	NTL_EXEC_RANGE_END;
	EvalMultMany_inpace_modImprove(perm_share_final, relin_keys, seal_context, bfv_secret_key);

	// perms[0] = EvalMultMany_inpace_modImprove(perms, relin_keys, seal_context, bfv_secret_key);

	time_end = chrono::high_resolution_clock::now();
	cout << "** [TIME] ** EvalMultMany_inpace_modImprove: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
	total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
	
	cout << "--[NOISE]-- after multiply to single perm vector: " << decryptor.invariant_noise_budget(perm_share_final[0]) << endl;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	time_start = chrono::high_resolution_clock::now();
    Ciphertext final_perm_vec = perm_share_final[0];
    // print_ct_to_pl(final_perm_vec, seal_context, bfv_secret_key, ring_dim);
    
    // oblivious expand the result into ring_dim ciphertexts, each encode the 0 or token value as the constant
    vector<Ciphertext> expanded = expand(seal_context, bfv_params, bfv_secret_key, final_perm_vec, ring_dim, gal_keys_expand);
    // print_ct_to_pl(expanded[2], seal_context, bfv_secret_key, ring_dim);

	time_end = chrono::high_resolution_clock::now();
	cout << "** [TIME] ** expand: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
	total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // multiply the plaintext token together with the expanded 0/1 plaintext vector
	time_start = chrono::high_resolution_clock::now();
    Ciphertext result;

    for (int i = 0; i < 4; i++) {
		for (int i = 0; i < ring_dim; i++) {
			if (i == 0) {
				evaluator.multiply_plain(expanded[i], tokens[i], result);
			} else {
				Ciphertext tmp;
				evaluator.multiply_plain(expanded[i], tokens[i], tmp);
				evaluator.add_inplace(result, tmp);
			}
		}
    }
    
    time_end = chrono::high_resolution_clock::now();
	cout << "** [TIME] ** multiply with pk: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
	cout << "****** TOTAL time: " << total_time << endl;


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    cout << "--[NOISE]-- final noise: " << decryptor.invariant_noise_budget(result) << endl;
    // print_ct_to_pl(result, seal_context, bfv_secret_key, ring_dim);


















    // // multiply all perm vec ciphertexts, log(n) depth of ct-multi
    // EvalMultMany_inpace_modImprove(perms, relin_keys, seal_context, bfv_secret_key);

    // Ciphertext final_perm_vec = perms[0];
    // print_ct_to_pl(final_perm_vec, seal_context, bfv_secret_key, ring_dim);
    // cout << "--[NOISE]-- after multiply to single perm vector: " << decryptor.invariant_noise_budget(final_perm_vec) << endl;


    // // inverse_slotToCoeff perm vector into slot values, preparing for multiplying with tokens
    // int sq_ct = sqrt(ring_dim/2);
    // vector<Ciphertext> packSIC_sqrt_list(2*sq_ct);
    // Ciphertext final_copy(final_perm_vec);
    // evaluator.rotate_columns_inplace(final_copy, gal_keys_slotToCoeff);

    // packSIC_sqrt_list[0] = final_perm_vec;
    // packSIC_sqrt_list[sq_ct] = final_copy;

    // for (int c = 1; c < sq_ct; c++) {
    //   evaluator.rotate_rows(packSIC_sqrt_list[c-1], sq_ct, gal_keys_slotToCoeff, packSIC_sqrt_list[c]);
    //   evaluator.rotate_rows(packSIC_sqrt_list[c-1+sq_ct], sq_ct, gal_keys_slotToCoeff, packSIC_sqrt_list[c+sq_ct]);
    // }
    // for (int c = 0; c < sq_ct; c++) {
    //   evaluator.transform_to_ntt_inplace(packSIC_sqrt_list[c]);
    //   evaluator.transform_to_ntt_inplace(packSIC_sqrt_list[c+sq_ct]);
    // }
    // Ciphertext perm_val = inverse_slotToCoeff_WOPrepreocess(seal_context, seal_context, packSIC_sqrt_list,
    // 							    gal_keys_slotToCoeff, ring_dim, p);
    // cout << "--[NOISE]-- after coeff to slot: " << decryptor.invariant_noise_budget(perm_val) << endl;

    // print_ct_to_vec(perm_val, seal_context, bfv_secret_key, ring_dim);


    // // multiply with the tokens (represented also in slot value)
    // vector<uint64_t> token(ring_dim);
    // for (int i = 0; i < (int) token.size(); i++) {
    //   token[i] = random_uint64() % 65537;
    // }
    // cout << "==================================== Token ====================================\n";
    // for (int i = 0; i < group_size; i++) {
    //   for (int j = 0; j < batch_size; j++) {
    // 	cout << token[j*group_size + i] << " ";
    //   }
    //   cout << endl;
    // }
    // cout << "===============================================================================\n";
    // Plaintext token_pl;
    // batch_encoder.encode(token, token_pl);
    // evaluator.multiply_plain_inplace(perm_val, token_pl);
    // // print_ct_to_vec(perm_val, seal_context, bfv_secret_key, ring_dim);
    // cout << "--[NOISE]-- after multiply plain with token: " << decryptor.invariant_noise_budget(perm_val) << endl;

    // // slotToCoeff, prepare for oblivious expansion which happens on the plaintext level
    // uint64_t inv = modInverse(ring_dim, p);
    // Ciphertext perm_copy(perm_val);
    // evaluator.rotate_columns_inplace(perm_copy, gal_keys_slotToCoeff);

    // packSIC_sqrt_list[0] = perm_val;
    // packSIC_sqrt_list[sq_ct] = perm_copy;

    // for (int c = 1; c < sq_ct; c++) {
    //   evaluator.rotate_rows(packSIC_sqrt_list[c-1], sq_ct, gal_keys_slotToCoeff, packSIC_sqrt_list[c]);
    //   evaluator.rotate_rows(packSIC_sqrt_list[c-1+sq_ct], sq_ct, gal_keys_slotToCoeff, packSIC_sqrt_list[c+sq_ct]);
    // }
    // for (int c = 0; c < sq_ct; c++) {
    //   evaluator.transform_to_ntt_inplace(packSIC_sqrt_list[c]);
    //   evaluator.transform_to_ntt_inplace(packSIC_sqrt_list[c+sq_ct]);
    // }
    // Ciphertext selected_token = slotToCoeff_WOPrepreocess(seal_context, seal_context, packSIC_sqrt_list,
    // 							  gal_keys_slotToCoeff, ring_dim, p, inv);

    // // oblivious expand the result into ring_dim ciphertexts, each encode the 0 or token value as the constant
    // vector<Ciphertext> expanded = expand(seal_context, bfv_params, bfv_secret_key, selected_token, ring_dim, gal_keys_expand);
    // print_ct_to_pl(expanded[2], seal_context, bfv_secret_key, ring_dim);
    // print_ct_to_pl(expanded[10], seal_context, bfv_secret_key, ring_dim);
    // print_ct_to_pl(expanded[18], seal_context, bfv_secret_key, ring_dim);
    // print_ct_to_pl(expanded[26], seal_context, bfv_secret_key, ring_dim);

    // vector<Ciphertext> rotated(batch_size);
    // rotated[0] = expanded[2];
    // for (int i = 0; i < ring_dim; i++) {
    //   pp.data()[i] = 0;
    // }
    // pp.data()[1] = 1;
    // evaluator.multiply_plain(expanded[10], pp, rotated[1]);

    // for (int i = 0; i < ring_dim; i++) {
    //   pp.data()[i] = 0;
    // }
    // pp.data()[2] = 1;
    // evaluator.multiply_plain(expanded[18], pp, rotated[2]);

    // for (int i = 0; i < ring_dim; i++) {
    //   pp.data()[i] = 0;
    // }
    // pp.data()[3] = 1;
    // evaluator.multiply_plain(expanded[26], pp, rotated[3]);

    // evaluator.add_inplace(rotated[0], rotated[1]);
    // evaluator.add_inplace(rotated[2], rotated[3]);
    // evaluator.add_inplace(rotated[0], rotated[2]);
    
    // // Ciphertext test;
    // // for (int i = 0; i < ring_dim; i++) {
    // //   pp.data()[i] = 0;
    // // }
    // // evaluator.add_plain(expanded[0], pp, test);
    // // // evaluator.add(expanded[1], expanded[2], test);

    // // // Plaintext pp;
    // // for (int i = 0; i < ring_dim; i++) {
    // //   pp.data()[i] = 0;
    // // }
    // // pp.data()[0] = 1;
    // // evaluator.multiply_plain_inplace(test, pp);
    // // print_ct_to_pl(test, seal_context, bfv_secret_key, ring_dim);

    // // // vector<Ciphertext> selected_token_batch(batch_size);
    // // // for (int i = 0; i < batch_size; i++) {
    // // //   selected_token_batch[i] = expanded[i * group_size];
    // // //   for (int j = 1; j < group_size; j++) {
    // // //   	evaluator.add_inplace(selected_token_batch[i], expanded[i * group_size + j]);
    // // //   }
    // // // }

    // cout << "--[NOISE]-- after oblivious expansion: " << decryptor.invariant_noise_budget(rotated[0]) << endl;

    // time_end = chrono::high_resolution_clock::now();
    // cout << "*** TOTAL time: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;

    // print_ct_to_pl(rotated[0], seal_context, bfv_secret_key, ring_dim);
    return 0;
}
