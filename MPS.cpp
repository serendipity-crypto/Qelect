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

	int numcores = 1;

    int group_size = 256;
    int batch_size = ring_dim / group_size;
    // int batch_size = ring_dim / group_size;

    EncryptionParameters bfv_params(scheme_type::bfv);
    bfv_params.set_poly_modulus_degree(ring_dim);
    auto coeff_modulus = CoeffModulus::Create(ring_dim, { 30, 60, 60, 60,  
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


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
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
    	perm_share_final[i] = EvalMultMany_inpace_modImprove(perms, relin_keys, seal_context, bfv_secret_key);
    }
	NTL_EXEC_RANGE_END;
	Ciphertext final_perm_vec = EvalMultMany_inpace_modImprove(perm_share_final, relin_keys, seal_context, bfv_secret_key);

	time_end = chrono::high_resolution_clock::now();
	cout << "** [TIME] ** EvalMultMany_inpace_modImprove: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
	total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
	
	cout << "--[NOISE]-- after multiply to single perm vector: " << decryptor.invariant_noise_budget(perm_share_final[0]) << endl;


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	time_start = chrono::high_resolution_clock::now();
    vector<Ciphertext> expanded_subtree_roots = subExpand(bfv_secret_key, seal_context, bfv_params, final_perm_vec,
                                                          ring_dim, gal_keys_expand, numcores);

    print_ct_to_pl(expanded_subtree_roots[0], seal_context, bfv_secret_key, ring_dim);
    
    // oblivious expand the result into ring_dim ciphertexts, each encode the 0 or token value as the constant

    vector<vector<Ciphertext>> expanded(numcores, vector<Ciphertext>(ring_dim/numcores));

    NTL_EXEC_RANGE(numcores, first, last);
    for (int i = first; i < last; i++) {
        expanded[i] = expand(seal_context, bfv_params, expanded_subtree_roots[i], ring_dim, gal_keys_expand,
                             ring_dim/numcores);
    }
    NTL_EXEC_RANGE_END;

    // print_ct_to_pl(expanded[0][0], seal_context, bfv_secret_key, ring_dim);


	time_end = chrono::high_resolution_clock::now();
	cout << "** [TIME] ** expand: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
	total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // multiply the plaintext token together with the expanded 0/1 plaintext vector
	time_start = chrono::high_resolution_clock::now();
    vector<Ciphertext> result_tmp(numcores);

    // cout << "--[NOISE]-- expanded noise: " << decryptor.invariant_noise_budget(expanded[0][0]) << endl;
    // cout << "--[NOISE]-- initial result noise: " << decryptor.invariant_noise_budget(result) << endl;

    NTL_EXEC_RANGE(numcores, first, last);
    int batch_share = ring_dim/numcores;
    for (int t = first; t < last; t++) {
        for (int i = t * batch_share; i < (t+1) * batch_share; i++) {
            if (i % batch_share == 0) { // skip the first one, done outside the loop already
                evaluator.multiply_plain(expanded[t][0], tokens[i], result_tmp[t]);
            } else {
                Ciphertext tmp;
                evaluator.multiply_plain(expanded[t][i % batch_share], tokens[i], tmp);
                evaluator.add_inplace(result_tmp[t], tmp);
            }
        }
    }
    NTL_EXEC_RANGE_END;

    // cout << "--[NOISE]-- final noise: " << decryptor.invariant_noise_budget(result_tmp[0]) << endl;

    for (int i = 1; i < (int) result_tmp.size(); i++) {
        evaluator.add_inplace(result_tmp[0], result_tmp[i]);
    }
    
    time_end = chrono::high_resolution_clock::now();
	cout << "** [TIME] ** multiply with pk: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
	cout << "****** TOTAL time: " << total_time << endl;


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    cout << "--[NOISE]-- final noise: " << decryptor.invariant_noise_budget(result_tmp[0]) << endl;
    // print_ct_to_pl(result, seal_context, bfv_secret_key, ring_dim);

    return 0;
}
