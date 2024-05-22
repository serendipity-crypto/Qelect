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
  
    int group_size = 8;
    
    int ring_dim = group_size; // for 200 people, can encrypt ring_dim / 200 Z_p element
    int n = 4;
    int p = 65537;

    int numcores = 1;

    EncryptionParameters bfv_params(scheme_type::bfv);
    bfv_params.set_poly_modulus_degree(ring_dim);
    auto coeff_modulus = CoeffModulus::Create(ring_dim, { 
                                                          60, 45, 60, 60,
                                                          60, 60, 60, 60,
                                                          60, 40, 60, 60
                                                        //   60, 60, 60, 60
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

    GaloisKeys gal_keys, gal_keys_expand, gal_keys_coeff;

    vector<int> stepsfirst = {1};
    keygen.create_galois_keys(stepsfirst, gal_keys);

    // gal keys for slotToCoeff
    vector<Modulus> coeff_modulus_coeff = coeff_modulus;
	coeff_modulus_coeff.erase(coeff_modulus_coeff.begin() + 10, coeff_modulus_coeff.end()-1);
	EncryptionParameters parms_coeff = bfv_params;
	parms_coeff.set_coeff_modulus(coeff_modulus_coeff);
	SEALContext context_coeff = SEALContext(parms_coeff, true, sec_level_type::none);
	Evaluator evaluator_next(context_coeff);

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
    
    Ciphertext perm;
    Plaintext pp;
    vector<uint64_t> perm_v(ring_dim);
    for (int g = 0; g < group_size; g++) {
        for (int i = 0; i < ring_dim; i++) {
            perm_v[i] = (i*49153) % p;
        }

        batch_encoder.encode(perm_v, pp);
        encryptor.encrypt(pp, perm);
        if (g == 0) cout << "-- [Noise] -- Initial noise: " << decryptor.invariant_noise_budget(perm) << endl;
        saveCiphertext(perm, g);
    }

    // prepare ring_dim different tokens, each is of size ring_dim, in plaintext form
    // for simplicity, use just one token
    Plaintext tokens;
    tokens.resize(ring_dim);
    tokens.parms_id() = parms_id_zero;
    for (int j = 0; j < (int) ring_dim; j++) {
        tokens.data()[j] = random_uint64() % 65537;
    }
    cout << "After param generation." << endl;

    chrono::high_resolution_clock::time_point time_start, time_end;
    uint64_t total_time = 0;
    time_start = chrono::high_resolution_clock::now();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // add all perm vec ciphertexts

    NTL::SetNumThreads(numcores);
    vector<Ciphertext> perm_share_final(numcores);
    uint64_t batch_perm_size = group_size / numcores;
    NTL_EXEC_RANGE(numcores, first, last);
    for (int i = first; i < last; i++) {
        perm_share_final[i] =  EvalAddMany_inpace_modImprove_extract_load(i*batch_perm_size, seal_context, 
                                                                          bfv_secret_key, batch_perm_size);
    }
    NTL_EXEC_RANGE_END;
    Ciphertext final_perm_vec = EvalAddMany_inpace_modImprove_extract_multi_core(perm_share_final, seal_context,
                                                                                 bfv_secret_key);

    cout << "-- [Noise] -- After adding up all perm: " << decryptor.invariant_noise_budget(final_perm_vec) << endl;

    time_end = chrono::high_resolution_clock::now();
    cout << "** [TIME] ** EvalMultMany_inpace_modImprove: " << \
          chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    evaluator.rotate_rows_inplace(final_perm_vec, 1, gal_keys); // place the second element to the first

    // extract the first element
    vector<uint64_t> extraction_v(ring_dim);
    for (int i = 0; i < ring_dim; i++) {
        extraction_v[i] = 0;
    }
    extraction_v[0] = 1;
    Plaintext extraction_p;
    batch_encoder.encode(extraction_v, extraction_p);
    evaluator.multiply_plain_inplace(final_perm_vec, extraction_p);

    evaluator.mod_switch_to_next_inplace(final_perm_vec);

    // slot_to_coeff so that the first slot is now serving as the constant term of the poly
    decryptor.decrypt(final_perm_vec, pp);
    batch_encoder.decode(pp, perm_v);
    cout << perm_v << endl;

    final_perm_vec = slotToCoeff_WOPrepreocess(context_coeff, final_perm_vec, gal_keys_coeff, ring_dim, p);

    decryptor.decrypt(final_perm_vec, pp);
    for (int i = 0; i < ring_dim; i++) {
        cout << pp.data()[i] << " ";
    }
    cout << endl;

    // compare with base_rot, and "blind rotate" it from "v" as a constant poly to "x^v" as a monomial
    vector<uint64_t> blind_rot_base_v(ring_dim);
    for (int i = 0; i < ring_dim; i++) {
        blind_rot_base_v[i] = p - i; // encrypt -i
    }
    Plaintext blind_rot_base_p;
    batch_encoder.encode(blind_rot_base_v, blind_rot_base_p);
    evaluator.add_plain_inplace(final_perm_vec, blind_rot_base_p);

    map<int, bool> raise_mod = {{4, false}, {16, false}, {64, false}, {256, false}};
    final_perm_vec = raisePowerToPrime(context_coeff, relin_keys_raise, final_perm_vec, raise_mod, raise_mod, 256, 256, p); // all one except v-th slot is zero


    vector<uint64_t> all_ones_v(ring_dim);
    for (int i = 0; i < ring_dim; i++) {
        all_ones_v[i] = 1; // encrypt 1
    }
    Plaintext all_ones_p;
    batch_encoder.encode(all_ones_v, all_ones_p);
    evaluator.negate_inplace(final_perm_vec);
    evaluator.add_plain_inplace(final_perm_vec, all_ones_p); // all zero except v-th slot is one

    // while(seal_context.last_parms_id() != final_perm_vec.parms_id()) {
	// 	evaluator.mod_switch_to_next_inplace(final_perm_vec);
	// }

    cout << "-- [Noise] -- After raise power for oblivious expansion: " << decryptor.invariant_noise_budget(final_perm_vec) << endl;

    decryptor.decrypt(final_perm_vec, pp);
    batch_encoder.decode(pp, perm_v);
    cout << perm_v << endl;

    final_perm_vec = slotToCoeff_WOPrepreocess(context_coeff, final_perm_vec, gal_keys_coeff, ring_dim, p);

    decryptor.decrypt(final_perm_vec, pp);
    for (int i = 0; i < ring_dim; i++) {
        cout << pp.data()[i] << " ";
    }
    cout << endl;


    vector<uint64_t> test_v(ring_dim);
    Plaintext test_p;
    Ciphertext test_c;
    for (int i = 0; i < ring_dim; i++) {
        test_v[i] = 0;
    }
    test_v[1] = 2;
    batch_encoder.encode(test_v, test_p);
    encryptor.encrypt(test_p, test_c);

    evaluator.mod_switch_to_next_inplace(test_c);

    decryptor.decrypt(test_c, pp);
    batch_encoder.decode(pp, perm_v);
    cout << perm_v << endl;

    test_c = slotToCoeff_WOPrepreocess(context_coeff, test_c, gal_keys_coeff, ring_dim, p);

    decryptor.decrypt(test_c, pp);
    for (int i = 0; i < ring_dim; i++) {
        cout << pp.data()[i] << " ";
    }
    cout << endl;


    // batch_encoder.decode(pp, perm_v);
    // cout << perm_v << endl;



    // perm_power = slotToCoeff_WOPrepreocess_time(seal_context, seal_context_next, perm_power,
    //                                             gal_keys_slotToCoeff, process_u_time[i], 128, ring_dim,
    //                                             t, inv);

    // vector<Ciphertext> expanded_subtree_roots_multi_core = subExpand(sk_expand, context_expand, parms_expand,
    //                                                              perm_power, ring_dim,
    //                                                              gal_keys_expand, numcores); 

    // // oblivious expand the result into ring_dim ciphertexts, each encode the 0 or token value as the constant
    // vector<vector<Ciphertext>> expanded_leaf(numcores, vector<Ciphertext>(ring_dim/numcores));

    // NTL_EXEC_RANGE(numcores, first, last);
    // for (int i = first; i < last; i++) {
    //     // expand each 1 out of the 8 subroots to leaf level
    //     expanded_leaf[i] = expand(context_expand, parms_expand, expanded_subtree_roots_multi_core[i], ring_dim,
    //                               gal_keys_expand, ring_dim/numcores);
    // }
    // NTL_EXEC_RANGE_END;
    // time_end = chrono::high_resolution_clock::now();
    // cout << "** [TIME] ** Second + Leaf level expand time: " << 
    //       chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    // total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();



    return 0;
}
