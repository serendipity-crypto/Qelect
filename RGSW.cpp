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
  
    int ring_dim = 2048; // for 200 people, can encrypt ring_dim / 200 Z_p element
    int n = 512;
    int p = 65537;

    EncryptionParameters bfv_params(scheme_type::bfv);
    bfv_params.set_poly_modulus_degree(ring_dim);
    auto coeff_modulus = CoeffModulus::Create(ring_dim, { 
                                                          60, 45, 60, 60,
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

    KeyGenerator keygen(seal_context);
    SecretKey bfv_secret_key = keygen.secret_key();
    PublicKey bfv_public_key;
    keygen.create_public_key(bfv_public_key);

    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(seal_context, bfv_public_key);
    Evaluator evaluator(seal_context);
    BatchEncoder batch_encoder(seal_context);
    Decryptor decryptor(seal_context, bfv_secret_key);


    vector<uint64_t> vv1(ring_dim), vv2(ring_dim);
    for (int i = 0; i < ring_dim; i++) {
        vv1[i] = 0;
        vv2[i] = 0;
    }
    vv1[10] = 1;
    vv2[1] = 1;


    Plaintext pp1, pp2, pp_sk;
    batch_encoder.encode(vv1, pp1);
    batch_encoder.encode(vv2, pp2);

    Ciphertext cc1, cc2_1, cc2_2, cc_sk;
    encryptor.encrypt(pp1, cc1); // m0
    encryptor.encrypt(pp2, cc2_1); // m1
    encryptor.encrypt(pp2, cc2_2); // m1


    pp_sk.resize(ring_dim);
    pp_sk.parms_id() = parms_id_zero;
    for (int i = 0; i < ring_dim; i++) { // prepare a ring elmt equals to 1
        if (i == 0) {
            pp_sk.data()[i] = 1;
        } else {
            pp_sk.data()[i] = 0;
        }
    }
    encryptor.encrypt_symmetric(pp_sk, cc_sk); // prepare encryption of -sk

    evaluator.multiply_inplace(cc2_1, cc_sk); // - sk*m
    evaluator.relinearize_inplace(cc2_1, relin_keys);

    evaluator.mod_switch_to_inplace(cc1, cc2_1.parms_id());
    evaluator.mod_switch_to_inplace(cc2_2, cc2_1.parms_id());

    vector<uint64_t> m0_a(ring_dim), m0_b(ring_dim), sm1_a(ring_dim), sm1_b(ring_dim), m1_a(ring_dim), m1_b(ring_dim);
    for (int i = 0; i < ring_dim; i++) {
        m0_a[i] = cc1.data(1)[i];
        m0_b[i] = cc1.data(0)[i];

        m1_a[i] = cc2_2.data(1)[i];
        m1_b[i] = cc2_2.data(0)[i];

        sm1_a[i] = cc2_1.data(1)[i];
        sm1_b[i] = cc2_1.data(0)[i];
    }


}