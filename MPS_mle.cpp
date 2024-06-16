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
  
    int ring_dim = 8; // for 200 people, can encrypt ring_dim / 200 Z_p element
    int n = 4;
    int p = 65537;

    // int numcores = 8;

    // int group_size = 4096;
    // int batch_size = ring_dim / sqrt(group_size);

    EncryptionParameters bfv_params(scheme_type::bfv);
    bfv_params.set_poly_modulus_degree(ring_dim);
    auto coeff_modulus = CoeffModulus::Create(ring_dim, { 
                                                          60,
                                                          60, 60, 60
                                                        });
    bfv_params.set_coeff_modulus(coeff_modulus);
    bfv_params.set_plain_modulus(p);

    cout << "after setting param.\n";

    prng_seed_type seed;
    for (auto &i : seed) {
        i = random_uint64();
    }
    auto rng = make_shared<Blake2xbPRNGFactory>(Blake2xbPRNGFactory(seed));
    bfv_params.set_random_generator(rng);

    SEALContext seal_context(bfv_params, true, sec_level_type::none);
    // primitive_root = seal_context.first_context_data()->plain_ntt_tables()->get_root();
    // cout << "primitive root: " << primitive_root << endl;

    KeyGenerator new_key_keygen(seal_context, n);
    SecretKey new_key = new_key_keygen.secret_key();
    inverse_ntt_negacyclic_harvey(new_key.data().data(), seal_context.key_context_data()->small_ntt_tables()[0]);

    KeyGenerator keygen(seal_context);
    SecretKey bfv_secret_key = keygen.secret_key();

    MemoryPoolHandle my_pool = MemoryPoolHandle::New();

    PublicKey bfv_public_key;
    keygen.create_public_key(bfv_public_key);

    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(seal_context, bfv_public_key);
    Evaluator evaluator(seal_context);
    BatchEncoder batch_encoder(seal_context);
    Decryptor decryptor(seal_context, bfv_secret_key);

    GaloisKeys gal_keys;

    // vector<Modulus> coeff_modulus_expand = coeff_modulus;
	// coeff_modulus_expand.erase(coeff_modulus_expand.begin() + 2, coeff_modulus_expand.end()-1);
	// EncryptionParameters parms_expand = bfv_params;
	// parms_expand.set_coeff_modulus(coeff_modulus_expand);
	// SEALContext context_expand = SEALContext(parms_expand, true, sec_level_type::none);

	// SecretKey sk_expand;
	// sk_expand.data().resize(coeff_modulus_expand.size() * ring_dim);
	// sk_expand.parms_id() = context_expand.key_parms_id();
	// util::set_poly(bfv_secret_key.data().data(), ring_dim, coeff_modulus_expand.size() - 1, sk_expand.data().data());
	// util::set_poly(
	// 		bfv_secret_key.data().data() + ring_dim * (coeff_modulus.size() - 1), ring_dim, 1,
	// 		sk_expand.data().data() + ring_dim * (coeff_modulus_expand.size() - 1));
	// KeyGenerator keygen_expand(context_expand, sk_expand);
	vector<uint32_t> galois_elts;
	for (int i = 0; i < ceil(log2(ring_dim)); i++) {
        galois_elts.push_back((ring_dim + exponentiate_uint(2, i)) / exponentiate_uint(2, i));
    }
	keygen.create_galois_keys(galois_elts, gal_keys);

    cout << "Gal key generated.\n";


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // prepare 256 permutation vectors, representing the first row of a permutation matrix
    Ciphertext perm;
    // vector<Ciphertext> perms(group_size);



    Plaintext pp;
    pp.resize(ring_dim);
    pp.parms_id() = parms_id_zero;
    // chrono::high_resolution_clock::time_point s1, e1;
    for (int i = 0; i < (int) ring_dim; i++) {
        pp.data()[i] = 0;
    }
    pp.data()[0] = 1;
    pp.data()[3] = 1;
    encryptor.encrypt(pp, perm);

    vector<Ciphertext> expanded_leaf(ring_dim);
    expanded_leaf = expand(seal_context, bfv_params, perm, ring_dim, gal_keys, ring_dim);

    for (int i = 0; i < (int) ring_dim; i++) {
        print_ct_to_pl(expanded_leaf[i], seal_context, bfv_secret_key, ring_dim);
    }


    return 0;
}
