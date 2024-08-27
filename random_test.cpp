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

    size_t poly_modulus_degree = 32768;
	int t = 65537;

	// step 3. generate detection key
	// recipient side
	EncryptionParameters parms(scheme_type::bfv);
	auto degree = poly_modulus_degree;
	parms.set_poly_modulus_degree(poly_modulus_degree);

	auto coeff_modulus = CoeffModulus::Create(poly_modulus_degree, { 60, 60,60, 60, 60});
	parms.set_coeff_modulus(coeff_modulus);
	parms.set_plain_modulus(t);

	prng_seed_type seed;
	for (auto &i : seed) {
		i = random_uint64();
	}
	auto rng = make_shared<Blake2xbPRNGFactory>(Blake2xbPRNGFactory(seed));
	parms.set_random_generator(rng);

	SEALContext context(parms, true, sec_level_type::none);
	KeyGenerator keygen(context);
	SecretKey secret_key = keygen.secret_key();

	PublicKey public_key;
	keygen.create_public_key(public_key);
	RelinKeys relin_keys;
	keygen.create_relin_keys(relin_keys);
	Encryptor encryptor(context, public_key);
	Evaluator evaluator(context);
	Decryptor decryptor(context, secret_key);
	BatchEncoder batch_encoder(context);


	GaloisKeys gal_keys, gal_keys_slotToCoeff, gal_keys_expand;
	vector<int> stepsfirst = {1};
	// only one rot key is needed for full level
	keygen.create_galois_keys(stepsfirst, gal_keys);

	/////////////////////////////////////////////////////////////// Rot Key gen
	vector<int> steps = {0};
	for(int i = 1; i < int(poly_modulus_degree/2); i *= 2) {
		steps.push_back(i);
	}


	//////////////////////////////////////////////////////
	vector<Modulus> coeff_modulus_expand = coeff_modulus;
	coeff_modulus_expand.erase(coeff_modulus_expand.begin() + 2, coeff_modulus_expand.end()-1);
	EncryptionParameters parms_expand = parms;
	parms_expand.set_coeff_modulus(coeff_modulus_expand);
	SEALContext context_expand = SEALContext(parms_expand, true, sec_level_type::none);

	SecretKey sk_expand;
	sk_expand.data().resize(coeff_modulus_expand.size() * degree);
	sk_expand.parms_id() = context_expand.key_parms_id();
	util::set_poly(secret_key.data().data(), degree, coeff_modulus_expand.size() - 1, sk_expand.data().data());
	util::set_poly(
			secret_key.data().data() + degree * (coeff_modulus.size() - 1), degree, 1,
			sk_expand.data().data() + degree * (coeff_modulus_expand.size() - 1));
	KeyGenerator keygen_expand(context_expand, sk_expand);
	vector<uint32_t> galois_elts;
	auto n = poly_modulus_degree;
	for (int i = 0; i < ceil(log2(poly_modulus_degree)); i++) {
		galois_elts.push_back((n + exponentiate_uint(2, i)) / exponentiate_uint(2, i));
	}
	keygen_expand.create_galois_keys(galois_elts, gal_keys_expand);

	// PublicKey public_key_last;
	// keygen_next.create_public_key(public_key_last);


	Plaintext pl;
	pl.resize(poly_modulus_degree);
	pl.parms_id() = parms_id_zero;
	for (int i = 0; i < (int) poly_modulus_degree; i++) {
		pl.data()[i] = 0;
	}
	pl.data()[2] = modInverse(32768, t);
	Ciphertext ran;
	encryptor.encrypt(pl, ran);

	evaluator.mod_switch_to_next_inplace(ran);
    evaluator.mod_switch_to_next_inplace(ran);

	cout << ran.coeff_modulus_size() << endl;

	vector<Ciphertext> expanded_subtree_leaves = subExpand(sk_expand, context_expand, parms_expand, ran, poly_modulus_degree, gal_keys_expand, 8);

	chrono::high_resolution_clock::time_point time_start, time_end;

	time_start = chrono::high_resolution_clock::now();
	
	// uint64_t inv = modInverse(32768, 65537);

	NTL::SetNumThreads(8);
	vector<vector<Ciphertext>> leaf_core(8, vector<Ciphertext>(32768/8));
	NTL_EXEC_RANGE(8, first, last);
    for (int i = first; i < last; i++) {
		leaf_core[i] = expand(context_expand, parms_expand, expanded_subtree_leaves[i], poly_modulus_degree, gal_keys_expand, poly_modulus_degree/8);
		

		for (int j = 0; j < (int) leaf_core[i].size(); j++) {
			Plaintext tmp;
			tmp.resize(poly_modulus_degree);
			tmp.parms_id() = parms_id_zero;
			for (int k = 0; k < (int) poly_modulus_degree; k++) {
				tmp.data()[k] = 0;
			}
			tmp.data()[0] = (i+1) % t;
			// tmp.data()[0] = 0;
            evaluator.multiply_plain_inplace(leaf_core[i][j], tmp);

			// evaluator.multiply_plain_inplace(leaf_core[i][j], pl);
        }

	}
	NTL_EXEC_RANGE_END;

	for (int iter = 0; iter < 8; iter++) {
		for (int i = 1; i < (int) leaf_core[0].size(); i++) {
			evaluator.add_inplace(leaf_core[iter][0], leaf_core[iter][i]);
		}
		if (iter !=0) {
			evaluator.add_inplace(leaf_core[0][0], leaf_core[iter][0]);
		}
	}

	cout << decryptor.invariant_noise_budget(leaf_core[0][0]) << endl;

	print_ct_to_vec(leaf_core[0][0], context_expand, sk_expand, 100);

	time_end = chrono::high_resolution_clock::now();
	cout << "time: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;

    // while(seal_context.last_parms_id() != test_ct.parms_id()) {
    //         evaluator.mod_switch_to_next_inplace(test_ct);
    // }

    // stringstream ss;
    // test_ct.save(ss);
    // ofstream datafile;
    // datafile.open ("../data/random_ct.txt");
    // datafile << ss.rdbuf();
    // datafile.close();



    // ifstream datafile;
    // datafile.open ("../data/random_ct.txt");

    // stringstream buf;
    // buf << datafile.rdbuf();

    // test_ct.load(seal_context, buf);

    
    // decryptor.decrypt(test_ct, pl);
    // batch_encoder.decode(pl, vvv);
    // cout << vvv << endl;

    // // for (int i = 0; i < ring_dim; i++) {
    // //     test_ct.data(0)[i] = random_uint64() % small_p;
    // // }
    // // for (int i = 0; i < ring_dim; i++) {
    // //     test_ct.data(1)[i] = random_uint64() % small_p;
    // // }

    // cout << decryptor.invariant_noise_budget(test_ct) << endl;



    // vector<uint64_t> extraction_v(ring_dim);
    // for (int i = 0; i < ring_dim; i++) {
    //     extraction_v[i] = 0;
    // }
    // extraction_v[0] = 1;
    // Plaintext extraction_p;
    // batch_encoder.encode(extraction_v, extraction_p);
    // evaluator.multiply_plain_inplace(test_ct, extraction_p);


    // decryptor.decrypt(test_ct, pl);
    // batch_encoder.decode(pl, vvv);
    // cout << vvv << endl;

    // decryptor.decrypt(test_ct, pl);
    // batch_encoder.decode(pl, vvv);
    // cout << vvv << endl;

    

    // evaluator.mod_switch_to_next_inplace(test_ct);

    // cout << decryptor.invariant_noise_budget(test_ct) << endl;

    // decryptor.decrypt(test_ct, pl);
    // batch_encoder.decode(pl, vvv);
    // cout << vvv << endl;


    // pl.resize(ring_dim);
    // pl.parms_id() = parms_id_zero;
    // for (int i = 0; i < ring_dim; i++) {
    //     pl.data()[i] = 10*i;
    // }
    // Ciphertext test_ct;
    // encryptor.encrypt(pl, test_ct);

    // while(seal_context.last_parms_id() != test_ct.parms_id()) {
    //         evaluator.mod_switch_to_next_inplace(test_ct);
    // }

    // cout << endl << "b part:\n";
    // for (int i = 0; i < ring_dim; i++) {
    //     cout << test_ct.data(0)[i] << " ";
    // }
    // cout << endl << "a part:\n";
    // for (int i = 0; i < ring_dim; i++) {
    //     cout << test_ct.data(1)[i] << " ";
    // }

    // cout << endl << "sk:\n";
    // inverse_ntt_negacyclic_harvey(bfv_secret_key.data().data(), seal_context.key_context_data()->small_ntt_tables()[0]);
	// for (int i = 0; i < ring_dim; i++) {
	// 	cout << bfv_secret_key.data()[i] << " ";
	// }
	// cout << endl;


}