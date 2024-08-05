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

    EncryptionParameters bfv_params(scheme_type::bfv);
    bfv_params.set_poly_modulus_degree(ring_dim);
    auto coeff_modulus = CoeffModulus::Create(ring_dim, { 
                                                          60, 60, 60, 60, 60, 60,
                                                          60, 60, 60, 60, 60, 60,
                                                          60, 60, 60, 60,
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
        stepsfirst.push_back(- (1<<i));
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
    // coeff_modulus_expand.erase(coeff_modulus_expand.begin() + 2, coeff_modulus_expand.end()-1);
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
    keygen.create_galois_keys(galois_elts, gal_keys_expand);


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

    // into [committee_size] different degree-one poly
    time_start = chrono::high_resolution_clock::now();
   

	Ciphertext tmp_random;
    vector<uint64_t> vvv(ring_dim);
    for (int i = 0; i < ring_dim; i++) {
		vvv[i] = random_uint64();
    }
	batch_encoder.encode(vvv, pp);
	encryptor.encrypt(pp, tmp_random);


	Ciphertext extracted;
   	for (int i = 0; i < ring_dim; i++) {
		vvv[i] = 0;
    }
	vvv[0] = 1; // extractor
	batch_encoder.encode(vvv, pp);
	evaluator.multiply_plain(tmp_random, pp, extracted);
    
    cout << decryptor.invariant_noise_budget(extracted) << endl;

	// expand random values to all slots
    expandFirstToAll(bfv_secret_key, seal_context, extracted, gal_keys, ring_dim);
    Ciphertext extracted_all_slots = extracted;

	// Ciphertext extracted_all_slots = slotToCoeff_WOPrepreocess(seal_context, extracted, gal_keys_coeff, ring_dim, p, scalar, numcores);
	time_end = chrono::high_resolution_clock::now();
	cout << "	first : slotToCoeff_WOPrepreocess: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
	cout << decryptor.invariant_noise_budget(extracted_all_slots) << endl;
	
	evaluator.negate_inplace(extracted_all_slots); // [-r, -r, ..., -r]
    map<int, bool> raise_mod_1 = {{2, false}, {8, false}, {32, false}, {128, false}};
	map<int, bool> raise_mod_2 = {{4, false}, {16, false}, {64, false}, {256, false}};

	Ciphertext evaluated1, evaluated2;
	for (int i = 0; i < ring_dim; i++) {
		vvv[i] = i;
    }
	batch_encoder.encode(vvv, pp);
	evaluator.add_plain(extracted_all_slots, pp, evaluated1); // [-r, -r, ..., -r] + [0, 1, ..., n-1]
	evaluated1 = raisePowerToPrime(seal_context, relin_keys, evaluated1, raise_mod_2, raise_mod_2, 256, 256, p);
	time_end = chrono::high_resolution_clock::now();
	cout << "	raisePower: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
	cout << decryptor.invariant_noise_budget(evaluated1) << endl;

	for (int i = 0; i < ring_dim; i++) {
		vvv[i] = i+ring_dim;
    }
	batch_encoder.encode(vvv, pp);
	evaluator.add_plain(extracted_all_slots, pp, evaluated2); // [-r, -r, ..., -r] + [0, 1, ..., n-1]
	evaluated2 = raisePowerToPrime(seal_context, relin_keys, evaluated2, raise_mod_2, raise_mod_2, 256, 256, p);
	time_end = chrono::high_resolution_clock::now();
	cout << "	raisePower: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
	cout << decryptor.invariant_noise_budget(evaluated2) << endl;

	evaluator.add_inplace(evaluated1, evaluated2); // binary with only ONE "1"
	evaluated1 = slotToCoeff_WOPrepreocess(seal_context, evaluated1, gal_keys_coeff, ring_dim, p, scalar, numcores);
	time_end = chrono::high_resolution_clock::now();
	cout << "	second: slotToCoeff_WOPrepreocess: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
	cout << decryptor.invariant_noise_budget(evaluated1) << endl;

	// evaluator.mod_switch_to_next_inplace(evaluated1);
	cout << evaluated1.coeff_modulus_size() << endl;
	cout << decryptor.invariant_noise_budget(evaluated1) << endl;

	vector<Ciphertext> expanded_leaf(ring_dim);
    expanded_leaf = expand(seal_context, bfv_params, evaluated1, ring_dim, gal_keys_expand, ring_dim);
	time_end = chrono::high_resolution_clock::now();
	cout << "	expand: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
	cout << decryptor.invariant_noise_budget(expanded_leaf[0]) << endl;

	Ciphertext result;
	for (int i = 0; i < ring_dim; i++) {
		if (i == 0) {
			evaluator.multiply_plain(expanded_leaf[i], pp, result);
		} else {
			evaluator.multiply_plain_inplace(expanded_leaf[i], pp);
			evaluator.add_inplace(result, expanded_leaf[i]);
		}
	}
	cout << decryptor.invariant_noise_budget(result) << endl;
	time_end = chrono::high_resolution_clock::now();

	cout << "Total time: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;




    return 0;
}
