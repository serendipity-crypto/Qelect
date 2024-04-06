#pragma once

#include "math/ternaryuniformgenerator.h"
#include "math/discreteuniformgenerator.h"
#include "math/discretegaussiangenerator.h"
#include "seal/seal.h"
#include <iostream>
#include <random>

using namespace std;
using namespace lbcrypto;
using namespace seal;


struct tFHEParam{
    int n;
    int p;
    double std_dev; // for each party
    int m;
    int t; // number of parties, we assume t-t tFHE
    int ring_dim;
    int batch;
    tFHEParam(){
        n = 512;
        p = 65537;
        std_dev = 1.3;
        m = 16000; 
        t = 256;
        ring_dim = 2048;
        batch = ring_dim / n;
    }
    tFHEParam(int n, int p, double std_dev, int m, int t, int ring_dim)
    : n(n), p(p), std_dev(std_dev), m(m), t(t), ring_dim(ring_dim)
    {}
};

typedef SEALContext Context;

typedef SecretKey SecretKeyShare;
typedef PublicKey PublicKeyShare;

Context setup_tfhe(tFHEParam& param) {
    EncryptionParameters bfv_params(scheme_type::bfv);
    bfv_params.set_poly_modulus_degree(param.ring_dim);
    auto coeff_modulus = CoeffModulus::Create(param.ring_dim, { 30, 60, 60, 60,  
                                                                60, 60, 60, 60,
                                                                60 });
    bfv_params.set_coeff_modulus(coeff_modulus);
    bfv_params.set_plain_modulus(param.p);

    prng_seed_type seed;
    for (auto &i : seed) {
        i = random_uint64();
    }
    auto rng = make_shared<Blake2xbPRNGFactory>(Blake2xbPRNGFactory(seed));
    bfv_params.set_random_generator(rng);

    Context context(bfv_params, true, sec_level_type::none);
    return context;
}


SecretKeyShare generate_tfhe_sk_share(KeyGenerator& key_gen) {
    return key_gen.secret_key();
}

PublicKeyShare generate_tfhe_pk_share(KeyGenerator& key_gen) {
    PublicKeyShare public_key_share;
    key_gen.create_public_key(public_key_share);

    return public_key_share;
}

PublicKey generate_tfhe_pk(vector<PublicKeyShare>& pk_shares, tFHEParam& param,
                           Context& context) {

    // generate a fresh new key
    KeyGenerator keygen(context, param.n);
    SecretKey dummy = keygen.secret_key();
    PublicKey aggregated_public_key;
    keygen.create_public_key(aggregated_public_key);

    uint64_t large_p = 0; // extract the key level prime
    inverse_ntt_negacyclic_harvey(dummy.data().data(), context.key_context_data()->small_ntt_tables()[0]);
    for (int i = 0; i < param.n; i++) {
        if (dummy.data()[i] > 1) { //p-1, which is -1, where p is prime of the last level in coeff modulus
            large_p = (uint64_t) dummy.data()[i];
            break;
        }
    }

    inverse_ntt_negacyclic_harvey(aggregated_public_key.data().data(),
                                  context.key_context_data()->small_ntt_tables()[0]);
    for (int i = 0; i < param.n; i++) {
        aggregated_public_key.data()[i] = 0; // erase all values
    }

    for (int i = 0; i < param.t; i++) {
        inverse_ntt_negacyclic_harvey(pk_shares[i].data().data(),
                                      context.key_context_data()->small_ntt_tables()[0]);
        for (int j = 0; j < param.n; j++) {
            aggregated_public_key.data()[j] = (uint64_t) (pk_shares[i].data()[j] + aggregated_public_key.data()[j]) % large_p;
        }
        seal::util::RNSIter key_rns(pk_shares[i].data().data(), param.ring_dim);
        ntt_negacyclic_harvey(key_rns,
                              context.key_context_data()->parms().coeff_modulus().size(), 
                              context.key_context_data()->small_ntt_tables());
    }

    seal::util::RNSIter pk_rns(aggregated_public_key.data().data(), param.ring_dim);
    ntt_negacyclic_harvey(pk_rns,
                          context.key_context_data()->parms().coeff_modulus().size(),
                          context.key_context_data()->small_ntt_tables());

    
    return aggregated_public_key;
}

vector<uint64_t> partial_decrypt_tfhe(SecretKeyShare& sk, Ciphertext& ct, tFHEParam& param,
                                      Context& context) {


    inverse_ntt_negacyclic_harvey(sk.data().data(), context.key_context_data()->small_ntt_tables()[0]);
    uint64_t large_p = 0; // extract the key level prime
    for (int i = 0; i < param.n; i++) {
        if (sk.data()[i] > 1) { //p-1, which is -1, where p is prime of the last level in coeff modulus
            large_p = (uint64_t) sk.data()[i];
            break;
        }
    }

    vector<uint64_t> results(param.n, 0);

    for (int i = 0; i < param.n; i++) {
        results[i] = (uint64_t) ct.data(0)[i]; // assign c0 part
    }

    // now we calculate (c1 * sk)
    for (int cnt = 0; cnt < param.n; cnt++) {
        int ind = 0;
        uint64_t tmp = 0;
        for (int i = cnt; i >= 0 && ind < param.n; i--) {
            uint64_t ttmp = (uint64_t) ((uint128_t) ((uint64_t) ct.data(1)[i] * (uint64_t) sk.data()[ind]) % large_p);
            tmp = (tmp + ttmp) % large_p;
            ind++;
        }

        for (int i = param.ring_dim-1; i > param.ring_dim - param.n + cnt && ind < param.n; i--) {
            uint64_t ttmp = (uint64_t) ((uint128_t) ((uint64_t) ct.data(1)[i] * (uint64_t) sk.data()[ind]) % large_p);
            tmp = (tmp + ttmp) % large_p;
            ind++;
        }
        results[cnt] = (tmp + results[cnt]) % large_p;
    }

    return results; // no rounding here
}
