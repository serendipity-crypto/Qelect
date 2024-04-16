#include <NTL/BasicThreadPool.h>
#include <NTL/ZZ.h>
#include <thread>
#include "seal/util/polyarithsmallmod.h"
#include "seal/seal.h"

using namespace seal::util;
using namespace std;
using namespace seal;


void print_ct_to_pl(Ciphertext& ct, SEALContext& context, SecretKey& sk, uint64_t ring_dim = 4) {
  Decryptor decryptor(context, sk);
  Plaintext pp;
  decryptor.decrypt(ct, pp);
  for (int i = 0; i < (int) ring_dim; i++) {
    cout << pp.data()[i] << " ";
  }
  cout << endl;
}


void print_ct_to_vec(Ciphertext& ct, SEALContext& context, SecretKey& sk, uint64_t ring_dim = 4) {
  Decryptor decryptor(context, sk);
  BatchEncoder batch_encoder(context);
  Plaintext pp;
  vector<uint64_t> msg(ring_dim);
  decryptor.decrypt(ct, pp);
  batch_encoder.decode(pp, msg);

  cout << msg << endl;
}


Ciphertext extract_and_multiply_multi_core(const RelinKeys &relin_keys, const SEALContext& context,
                                     Ciphertext& ct1, Ciphertext& ct2, const int group_size = 4096,
				     const int ring_dim = 32768, const int numcores = 8) {
    Evaluator evaluator(context);
    NTL::SetNumThreads(numcores);

    int sq_group_size = sqrt(group_size); // extract each sq_group_size

    // vector<Ciphertext> add_tmp(sq_group_size);
    Ciphertext final;
    Plaintext pl;
    pl.resize(ring_dim);
    pl.parms_id() = parms_id_zero;
    for (int j = 0; j < (int) ring_dim; j++) { 
        if (j < sq_group_size) {
            pl.data()[j] = 1;
        } else {
            pl.data()[j] = 0;
        }
    }
    Ciphertext tmp1, tmp2;
    evaluator.multiply_plain(ct1, pl, tmp1);
    evaluator.multiply_plain(ct2, pl, tmp2);
    evaluator.multiply(tmp1, tmp2, final);


    int batch_size = sq_group_size / numcores;

    NTL_EXEC_RANGE(numcores, first, last);
    for (int c = first; c < last; c++) {
        for (int i = batch_size*c; i < batch_size*(c+1); i++) {

            if (i == 0) continue; // skip the initial setup

            // prepare the extraction plaintext
            Plaintext pl_core;
            pl_core.resize(ring_dim);
            pl_core.parms_id() = parms_id_zero;
            for (int j = 0; j < (int) ring_dim; j++) { 
                if (j < sq_group_size*(i+1) && j >= sq_group_size*i) {
                    pl_core.data()[j] = 1;
                } else {
                    pl_core.data()[j] = 0;
                }
            }

            Ciphertext tmp1_core, tmp2_core;
            evaluator.multiply_plain(ct1, pl_core, tmp1_core);
            evaluator.multiply_plain(ct2, pl_core, tmp2_core);

            Ciphertext tmp;
            evaluator.multiply(tmp1_core, tmp2_core, tmp);
            // evaluator.relinearize_inplace(tmp, relin_keys);
            evaluator.add_inplace(final, tmp);
        }
    }
    NTL_EXEC_RANGE_END;

    evaluator.relinearize_inplace(final, relin_keys);
    return final;
}

Ciphertext extract_and_multiply(const RelinKeys &relin_keys, const SEALContext& context, Ciphertext& ct1,
			  Ciphertext& ct2, const int group_size = 4096, const int ring_dim = 32768) {
    Evaluator evaluator(context);

    int sq_group_size = sqrt(group_size); // extract each sq_group_size

    Plaintext pl;
    pl.resize(ring_dim);
    pl.parms_id() = parms_id_zero;
    // vector<Ciphertext> add_tmp(sq_group_size);
    Ciphertext final;

    for (int i = 0; i < sq_group_size; i++) {
        // prepare the extraction plaintext
        for (int j = 0; j < (int) ring_dim; i++) { 
            if (j < sq_group_size*(i+1) && j >= sq_group_size*i) {
                pl.data()[i] = 1;
            } else {
                pl.data()[i] = 0;
            }
        }

        Ciphertext tmp1, tmp2;
        evaluator.multiply_plain(ct1, pl, tmp1);
        evaluator.multiply_plain(ct2, pl, tmp2);

        if (i == 0){
            evaluator.multiply(tmp1, tmp2, final);
            // evaluator.relinearize_inplace(final, relin_keys);
        } else {
            Ciphertext tmp;
            evaluator.multiply(tmp1, tmp2, tmp);
            // evaluator.relinearize_inplace(tmp, relin_keys);
            evaluator.add_inplace(final, tmp);
        } 
    }

    evaluator.relinearize_inplace(final, relin_keys);
    return final;
}


inline
Ciphertext EvalMultMany_inpace_modImprove_extract(vector<Ciphertext>& ciphertexts, const RelinKeys &relin_keys,
						  const SEALContext& context, SecretKey& sk) {
    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);
    int counter = 0;

    while(ciphertexts.size() != 1){
        for(size_t i = 0; i < ciphertexts.size()/2; i++){
	    ciphertexts[i] = extract_and_multiply(relin_keys, context, ciphertexts[i], ciphertexts[ciphertexts.size()/2+i]);
            // evaluator.multiply_inplace(ciphertexts[i], ciphertexts[ciphertexts.size()/2+i]);
            // evaluator.relinearize_inplace(ciphertexts[i], relin_keys);
            if(counter & 1) {
                evaluator.mod_switch_to_next_inplace(ciphertexts[i]);
            }
        }
        if(ciphertexts.size()%2 == 0)
            ciphertexts.resize(ciphertexts.size()/2);
        else{ // if odd, take the last one and mod down to make them compatible                                                                                                                                        
            ciphertexts[ciphertexts.size()/2] = ciphertexts[ciphertexts.size()-1];
            if(counter & 1) {
                evaluator.mod_switch_to_next_inplace(ciphertexts[ciphertexts.size()/2]);
            }
            ciphertexts.resize(ciphertexts.size()/2+1);
        }
        counter += 1;
    }

    Ciphertext res = ciphertexts[0];
    return res;
}

inline
Ciphertext EvalMultMany_inpace_modImprove(vector<Ciphertext>& ciphertexts, const RelinKeys &relin_keys,
				    const SEALContext& context, SecretKey& sk) {
    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);
    int counter = 0;

    while(ciphertexts.size() != 1){
        for(size_t i = 0; i < ciphertexts.size()/2; i++){
            evaluator.multiply_inplace(ciphertexts[i], ciphertexts[ciphertexts.size()/2+i]);
            evaluator.relinearize_inplace(ciphertexts[i], relin_keys);
            if(counter & 1) {
                evaluator.mod_switch_to_next_inplace(ciphertexts[i]);
            }
        }
        if(ciphertexts.size()%2 == 0)
            ciphertexts.resize(ciphertexts.size()/2);
        else{ // if odd, take the last one and mod down to make them compatible
            ciphertexts[ciphertexts.size()/2] = ciphertexts[ciphertexts.size()-1];
            if(counter & 1) {
                evaluator.mod_switch_to_next_inplace(ciphertexts[ciphertexts.size()/2]);
            }
            ciphertexts.resize(ciphertexts.size()/2+1);
        }
        counter += 1;
    }

    Ciphertext res = ciphertexts[0];
    return res;
}


inline
long power(long x, long y, long m)
{
    if (y == 0)
        return 1;
    long p = power(x, y / 2, m) % m;
    p = (p * p) % m;
 
    return (y % 2 == 0) ? p : (x * p) % m;
}

inline
long modInverse(long a, long m)
{
    return power(a, m - 2, m);
}


vector<regevCiphertext> extractRLWECiphertextToLWECiphertext(Ciphertext& rlwe_ct, const int ring_dim = poly_modulus_degree_glb,
                                                             const int n = 1024, const int p = prime_p,
                                                             const uint64_t big_prime = 1152921504589938689) {
    vector<regevCiphertext> results(ring_dim);

    prng_seed_type seed;
    for (auto &i : seed) {
        i = random_uint64();
    }
    auto rng = std::make_shared<Blake2xbPRNGFactory>(Blake2xbPRNGFactory(seed));
    RandomToStandardAdapter engine(rng->create());
    uniform_int_distribution<uint32_t> dist(0, 100);

    for (int cnt = 0; cnt < ring_dim; cnt++) {
        results[cnt].a = NativeVector(n);
        int ind = 0;
        for (int i = cnt; i >= 0 && ind < n; i--) {
            float temp_f = ((float) rlwe_ct.data(1)[i]) * ((float) p) / ((long double) big_prime);
            uint32_t decimal = (temp_f - ((int) temp_f)) * 100;
            float rounding = dist(engine) < decimal ? 1 : 0;

            long temp = ((int) (temp_f + rounding)) % p;
            results[cnt].a[ind] = temp < 0 ? p + temp : temp;

            ind++;
        }

        for (int i = ring_dim-1; i > ring_dim - n + cnt && ind < n; i--) {
            float temp_f = ((float) rlwe_ct.data(1)[i]) * ((float) p) / ((long double) big_prime);
            uint32_t decimal = (temp_f - ((int) temp_f)) * 100;
            float rounding = dist(engine) < decimal ? 1 : 0;

            long temp = ((int) (temp_f + rounding)) % p;
            results[cnt].a[ind] = -temp < 0 ? p-temp : -temp;

            ind++;
        }

        float temp_f = ((float) rlwe_ct.data(0)[cnt]) * ((float) p) / ((long double) big_prime);
        uint32_t decimal = temp_f - ((int) temp_f) * 100;
        float rounding = dist(engine) < decimal ? 1 : 0;

        long temp = ((int) (temp_f + rounding)) % p;
        results[cnt].b = temp % ((int) p);
    }

    return results;
}



inline void multiply_power_of_X(EncryptionParameters& enc_param, const Ciphertext &encrypted, Ciphertext &destination,
                                uint32_t index) {

    auto coeff_mod_count = enc_param.coeff_modulus().size() - 1;
    auto coeff_count = enc_param.poly_modulus_degree();
    auto encrypted_count = encrypted.size();

    destination = encrypted;

    for (int i = 0; i < (int) encrypted_count; i++) {
        for (int j = 0; j < (int) coeff_mod_count; j++) {
            negacyclic_shift_poly_coeffmod(encrypted.data(i) + (j * coeff_count),
                                           coeff_count, index,
                                           enc_param.coeff_modulus()[j],
                                           destination.data(i) + (j * coeff_count));
        }
    }
}

inline vector<Ciphertext> subExpand(const SecretKey& sk, const SEALContext& context, EncryptionParameters& enc_param,
                                    const Ciphertext &encrypted, uint32_t m, const GaloisKeys& galkey,
                                    int first_expansion_size, int t = 65537) {

    Evaluator evaluator(context);
    Plaintext two("2");
    Decryptor decryptor(context, sk);

    int logFirst = ceil(log2(first_expansion_size));

    vector<int> galois_elts;

    for (int i = 0; i < ceil(log2(m)); i++) {
        galois_elts.push_back((m + exponentiate_uint(2, i)) / exponentiate_uint(2, i));
    }

    vector<Ciphertext> temp;
    temp.push_back(encrypted);
    Ciphertext tempctxt;
    Ciphertext tempctxt_rotated;
    Ciphertext tempctxt_shifted;
    Ciphertext tempctxt_rotatedshifted;

    for (int i = 0; i < logFirst; i++) {
        vector<Ciphertext> newtemp(temp.size() << 1);
        int index_raw = (m << 1) - (1 << i);
        int index = (index_raw * galois_elts[i]) % (m << 1);

        for (uint32_t a = 0; a < temp.size(); a++) {

            evaluator.apply_galois(temp[a], galois_elts[i], galkey, tempctxt_rotated);
            evaluator.add(temp[a], tempctxt_rotated, newtemp[a]);
            multiply_power_of_X(enc_param, temp[a], tempctxt_shifted, index_raw);
            multiply_power_of_X(enc_param, tempctxt_rotated, tempctxt_rotatedshifted, index);
            evaluator.add(tempctxt_shifted, tempctxt_rotatedshifted, newtemp[a + temp.size()]);
        }

        temp = newtemp;
    }

    vector<Ciphertext>::const_iterator first = temp.begin();
    vector<Ciphertext>::const_iterator last = temp.begin() + first_expansion_size;
    vector<Ciphertext> newVec(first, last);

    return newVec;
}

inline vector<Ciphertext> expand(const SEALContext& context, EncryptionParameters& enc_param, const Ciphertext &encrypted,
                                 uint32_t m, const GaloisKeys& galkey, int stepSize = 1, int t = 65537) {

    Evaluator evaluator(context);
    Plaintext two("2");

    int first_expansion_size = m / stepSize;
    int logFirst = ceil(log2(first_expansion_size));
    int logm = ceil(log2(m));

    vector<int> galois_elts;

    for (int i = 0; i < ceil(log2(m)); i++) {
        galois_elts.push_back((m + exponentiate_uint(2, i)) / exponentiate_uint(2, i));
    }

    vector<Ciphertext> temp;
    temp.push_back(encrypted);
    Ciphertext tempctxt;
    Ciphertext tempctxt_rotated;
    Ciphertext tempctxt_shifted;
    Ciphertext tempctxt_rotatedshifted;

    for (int i = logFirst; i < logm - 1; i++) {
        vector<Ciphertext> newtemp(temp.size() << 1);
        int index_raw = (m << 1) - (1 << i);
        int index = (index_raw * galois_elts[i]) % (m << 1);

        for (int a = 0; a < (int) temp.size(); a++) {

            evaluator.apply_galois(temp[a], galois_elts[i], galkey, tempctxt_rotated);

            evaluator.add(temp[a], tempctxt_rotated, newtemp[a]);
            multiply_power_of_X(enc_param, temp[a], tempctxt_shifted, index_raw);
            multiply_power_of_X(enc_param, tempctxt_rotated, tempctxt_rotatedshifted, index);

            evaluator.add(tempctxt_shifted, tempctxt_rotatedshifted, newtemp[a + temp.size()]);
        }

        temp = newtemp;
    }

    // Last step of the loop
    vector<Ciphertext> newtemp(temp.size() << 1);
    int index_raw = (m << 1) - (1 << (logm - 1));
    int index = (index_raw * galois_elts[logm - 1]) % (m << 1);

    for (uint32_t a = 0; a < temp.size(); a++) {
        if (a >= (m - (1 << (logm - 1)))) { // corner case.
            evaluator.multiply_plain(temp[a], two, newtemp[a]); // plain multiplication by 2.
        } else {
            evaluator.apply_galois(temp[a], galois_elts[logm - 1], galkey, tempctxt_rotated);
            evaluator.add(temp[a], tempctxt_rotated, newtemp[a]);
            multiply_power_of_X(enc_param, temp[a], tempctxt_shifted, index_raw);
            multiply_power_of_X(enc_param, tempctxt_rotated, tempctxt_rotatedshifted, index);
            evaluator.add(tempctxt_shifted, tempctxt_rotatedshifted, newtemp[a + temp.size()]);
        }
    }

    vector<Ciphertext>::const_iterator first = newtemp.begin();
    vector<Ciphertext>::const_iterator last = newtemp.begin() + stepSize;
    vector<Ciphertext> newVec(first, last);

    return newVec;
}



// for a tree with m leaf node, m >> stepSize, we first expand it to a subtree with m / stepSize leaf node
// (i.e., this subtree is the top of the whole tree)
// and then for each leaf node in this subtree, expand it into a small subtree with stepSize leaf node
inline vector<Ciphertext> expand_standalone(const SEALContext& context, EncryptionParameters& enc_param,
				                 const SecretKey& sk, const Ciphertext &encrypted, uint32_t m,
                                 const GaloisKeys& galkey, int t = 65537) {

    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);
    Plaintext two("2");
    int logm = ceil(log2(m));

    vector<int> galois_elts;

    for (int i = 0; i < ceil(log2(m)); i++) {
        galois_elts.push_back((m + exponentiate_uint(2, i)) / exponentiate_uint(2, i));
    }

    vector<Ciphertext> temp;
    temp.push_back(encrypted);
    Ciphertext tempctxt;
    Ciphertext tempctxt_rotated;
    Ciphertext tempctxt_shifted;
    Ciphertext tempctxt_rotatedshifted;

    for (int i = 0; i < logm - 1; i++) {
        vector<Ciphertext> newtemp(temp.size() << 1);
        int index_raw = (m << 1) - (1 << i);
        int index = (index_raw * galois_elts[i]) % (m << 1);

        for (int a = 0; a < (int) temp.size(); a++) {
            evaluator.apply_galois(temp[a], galois_elts[i], galkey, tempctxt_rotated);

            evaluator.add(temp[a], tempctxt_rotated, newtemp[a]);
            multiply_power_of_X(enc_param, temp[a], tempctxt_shifted, index_raw);
            multiply_power_of_X(enc_param, tempctxt_rotated, tempctxt_rotatedshifted, index);

            evaluator.add(tempctxt_shifted, tempctxt_rotatedshifted, newtemp[a + temp.size()]);
        }

        temp = newtemp;
    }

    // Last step of the loop
    vector<Ciphertext> newtemp(temp.size() << 1);
    int index_raw = (m << 1) - (1 << (logm - 1));
    int index = (index_raw * galois_elts[logm - 1]) % (m << 1);

    for (uint32_t a = 0; a < temp.size(); a++) {
        if (a >= (m - (1 << (logm - 1)))) { // corner case.
            evaluator.multiply_plain(temp[a], two, newtemp[a]); // plain multiplication by 2.
        } else {
            evaluator.apply_galois(temp[a], galois_elts[logm - 1], galkey, tempctxt_rotated);
	    
            evaluator.add(temp[a], tempctxt_rotated, newtemp[a]);
            multiply_power_of_X(enc_param, temp[a], tempctxt_shifted, index_raw);
            multiply_power_of_X(enc_param, tempctxt_rotated, tempctxt_rotatedshifted, index);
            evaluator.add(tempctxt_shifted, tempctxt_rotatedshifted, newtemp[a + temp.size()]);
        }
    }

    vector<Ciphertext>::const_iterator first = newtemp.begin();
    vector<Ciphertext>::const_iterator last = newtemp.begin() + m;
    vector<Ciphertext> newVec(first, last);

    return newVec;
}

