
#include "seal/util/polyarithsmallmod.h"
#include "seal/seal.h"

using namespace seal::util;
using namespace std;
using namespace seal;


struct BootstrapParam {
    int ciphertextSpacePrime;
    int errorRange;
    int plaintextSpace;

    int firstLevelDegree; // used for CallToDegreeK in range check procedure
    int secondLevelDegree;

    int raisePower_firstLevel = 256;
    int raisePower_secondLevel = 256;

    BootstrapParam(int ciphertextSpacePrime, int errorRange, int plaintextSpace, int firstLevelDegree, int secondLevelDegree,
                   int raisePower_firstLevel = 256, int raisePower_secondLevel = 256)
    : ciphertextSpacePrime(ciphertextSpacePrime), errorRange(errorRange), plaintextSpace(plaintextSpace),
      firstLevelDegree(firstLevelDegree), secondLevelDegree(secondLevelDegree), raisePower_firstLevel(raisePower_firstLevel),
      raisePower_secondLevel(raisePower_secondLevel)
    {}
};


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

inline void calUptoDegreeK(vector<Ciphertext>& output, const Ciphertext& input, const int DegreeK, const RelinKeys &relin_keys,
                           const SEALContext& context, const bool skip_odd=false) {
    vector<int> calculated(DegreeK, 0);
    Evaluator evaluator(context);
    output[0] = input;
    calculated[0] = 1; // degree 1, x
    Ciphertext res, base;
    vector<int> numMod(DegreeK, 0);

    for(int i = DegreeK; i > 0; i--){
        if (skip_odd && i % 2 == 1) { // 0 is for degree 1, 1 is for degree 2, skip all 2k+1 degree
            calculated[i-1] = 1;
            output[i-1] = input;
        } else if(calculated[i-1] == 0){
            auto toCalculate = i;
            int resdeg = 0;
            int basedeg = 1;
            base = input;
            while(toCalculate > 0){
                if(toCalculate & 1){
                    toCalculate -= 1;
                    resdeg += basedeg;
                    if(calculated[resdeg-1] != 0){
                        res = output[resdeg - 1];
                    } else {
                        if(resdeg == basedeg){
                            res = base; // should've never be used as base should have made calculated[basedeg-1]
                        } else {
                            numMod[resdeg-1] = numMod[basedeg-1];

                            evaluator.mod_switch_to_inplace(res, base.parms_id()); // match modulus
                            evaluator.multiply_inplace(res, base);
                            evaluator.relinearize_inplace(res, relin_keys);
                            while(numMod[resdeg-1] < (ceil(log2(resdeg))/2)){
                                evaluator.mod_switch_to_next_inplace(res);
                                numMod[resdeg-1]+=1;
                            }
                        }
                        output[resdeg-1] = res;
                        calculated[resdeg-1] += 1;
                    }
                } else {
                    toCalculate /= 2;
                    basedeg *= 2;
                    if(calculated[basedeg-1] != 0){
                        base = output[basedeg - 1];
                    } else {
                        numMod[basedeg-1] = numMod[basedeg/2-1];
                        evaluator.square_inplace(base);
                        evaluator.relinearize_inplace(base, relin_keys);
                        while(numMod[basedeg-1] < (ceil(log2(basedeg))/2)){
                            evaluator.mod_switch_to_next_inplace(base);
                            numMod[basedeg-1]+=1;
                        }
                        output[basedeg-1] = base;
                        calculated[basedeg-1] += 1;
                    }
                }
            }
        }
    }

    for(size_t i = 0; i < output.size()-1; i++){
        evaluator.mod_switch_to_inplace(output[i], output[output.size()-1].parms_id()); // match modulus
    }
    return;
}


inline void calUptoDegreeK_bigPrime(vector<Ciphertext>& output, const Ciphertext& input, const int DegreeK, const RelinKeys &relin_keys,
                                    const SEALContext& context, map<int, bool>& modDownIndices, const bool skip_odd=false) {
    vector<int> calculated(DegreeK, 0);
    Evaluator evaluator(context);
    output[0] = input;
    calculated[0] = 1; // degree 1, x
    Ciphertext res, base;
    vector<int> numMod(DegreeK, 0);

    for(int i = DegreeK; i > 0; i--){
        if (skip_odd && i % 2 == 1) { // 0 is for degree 1, 1 is for degree 2, skip all 2k+1 degree
            calculated[i-1] = 1;
            output[i-1] = input;
        } else if(calculated[i-1] == 0){
            auto toCalculate = i;
            int resdeg = 0;
            int basedeg = 1;
            base = input;
            while(toCalculate > 0){
                if(toCalculate & 1){
                    toCalculate -= 1;
                    resdeg += basedeg;
                    if(calculated[resdeg-1] != 0){
                        res = output[resdeg - 1];
                    } else {
                        if(resdeg == basedeg){
                            res = base; // should've never be used as base should have made calculated[basedeg-1]
                        } else {
                            numMod[resdeg-1] = numMod[basedeg-1];

                            evaluator.mod_switch_to_inplace(res, base.parms_id()); // match modulus
                            evaluator.multiply_inplace(res, base);
                            evaluator.relinearize_inplace(res, relin_keys);
                            if(modDownIndices.count(resdeg) && !modDownIndices[resdeg]) {
                                modDownIndices[resdeg] = true;
                                evaluator.mod_switch_to_next_inplace(res);
                                numMod[resdeg-1]+=1;
                            }
                        }
                        output[resdeg-1] = res;
                        calculated[resdeg-1] += 1;
                    }
                } else {
                    toCalculate /= 2;
                    basedeg *= 2;
                    if(calculated[basedeg-1] != 0){
                        base = output[basedeg - 1];
                    } else {
                        numMod[basedeg-1] = numMod[basedeg/2-1];
                        evaluator.square_inplace(base);
                        evaluator.relinearize_inplace(base, relin_keys);
                        while(modDownIndices.count(basedeg) && !modDownIndices[basedeg]) {
                            modDownIndices[basedeg] = true;
                            evaluator.mod_switch_to_next_inplace(base);
                            numMod[basedeg-1]+=1;
                        }
                        output[basedeg-1] = base;
                        calculated[basedeg-1] += 1;
                    }
                }
            }
        }
    }

    for(size_t i = 0; i < output.size()-1; i++){
        evaluator.mod_switch_to_inplace(output[i], output[output.size()-1].parms_id()); // match modulus
    }
    return;
}



Ciphertext evaluateExtractedBFVCiphertext(const SEALContext& seal_context, vector<regevCiphertext>& lwe_ct_list,
                                          const vector<Ciphertext>& lwe_sk_sqrt_list, const GaloisKeys& gal_keys, const int lwe_sk_len,
                                          const vector<uint64_t>& q_shift_constant, const int degree = poly_modulus_degree_glb,
                                          const bool gateEval = false, const int q = prime_p) {
    Evaluator evaluator(seal_context);
    BatchEncoder batch_encoder(seal_context);

    // rotate sqrt(degree), and get sqrt(degree)'s lwe_sk_encrypted
    int sq_rt = sqrt(lwe_sk_len);
    vector<Ciphertext> result(sq_rt);
        
    for (int iter = 0; iter < sq_rt; iter++) {
        for (int j = 0; j < (int) lwe_sk_sqrt_list.size(); j++) {
            vector<uint64_t> lwe_ct_tmp(degree);
            for (int i = 0; i < degree; i++) {
                int ct_index = (i-iter) % (degree/2) < 0 ? (i-iter) % (degree/2) + degree/2 : (i-iter) % (degree/2);
                ct_index = i < degree/2 ? ct_index : ct_index + degree/2;
                int col_index = (i + j*sq_rt) % lwe_sk_len;
                lwe_ct_tmp[i] = lwe_ct_list[ct_index].a[col_index].ConvertToInt();
            }

            Plaintext lwe_ct_pl;
            batch_encoder.encode(lwe_ct_tmp, lwe_ct_pl);
            evaluator.transform_to_ntt_inplace(lwe_ct_pl, lwe_sk_sqrt_list[j].parms_id());

            if (j == 0) {
                evaluator.multiply_plain(lwe_sk_sqrt_list[j], lwe_ct_pl, result[iter]);
            } else {
                Ciphertext temp;
                evaluator.multiply_plain(lwe_sk_sqrt_list[j], lwe_ct_pl, temp);
                evaluator.add_inplace(result[iter], temp);
            }

        }
    }

    for (int i = 0; i < (int) result.size(); i++) {
        evaluator.transform_from_ntt_inplace(result[i]);
    }

    // sum up all sq_rt tmp results to the first one, each first rotate left one and add to the previous tmp result

    for (int iter = sq_rt-1; iter > 0; iter--) {
        evaluator.rotate_rows_inplace(result[iter], 1, gal_keys);
        evaluator.add_inplace(result[iter-1], result[iter]);
    }

    vector<uint64_t> b_parts(degree);
    for(int i = 0; i < degree; i++){
        b_parts[i] = lwe_ct_list[i].b.ConvertToInt();
    }

    Plaintext lwe_b_pl;
    batch_encoder.encode(b_parts, lwe_b_pl);
    evaluator.add_plain_inplace(result[0], lwe_b_pl);

    if (gateEval) {
        Plaintext q_shift_pl;
        batch_encoder.encode(q_shift_constant, q_shift_pl);
        evaluator.add_plain_inplace(result[0], q_shift_pl);
    }

    return result[0];
}


// assume lwe_sk_len is a power of 2, and has a square root
Ciphertext evaluatePackedLWECiphertext(const SEALContext& seal_context, vector<regevCiphertext>& lwe_ct_list,
                                       const vector<Ciphertext>& lwe_sk_sqrt_list, const GaloisKeys& gal_keys, const int lwe_sk_len,
                                       const vector<uint64_t>& q_shift_constant, const int degree = poly_modulus_degree_glb,
                                       const bool gateEval = false, const int q = prime_p) {
    Evaluator evaluator(seal_context);
    BatchEncoder batch_encoder(seal_context);

    // rotate sqrt(degree), and get sqrt(degree)'s lwe_sk_encrypted
    int sq_rt = sqrt(lwe_sk_len);
    vector<Ciphertext> result(sq_rt);
        
    for (int iter = 0; iter < sq_rt; iter++) {
        for (int j = 0; j < (int) lwe_sk_sqrt_list.size(); j++) {
            vector<uint64_t> lwe_ct_tmp(degree);
            for (int i = 0; i < degree; i++) {
                int ct_index = (i-iter) % (degree/2) < 0 ? (i-iter) % (degree/2) + degree/2 : (i-iter) % (degree/2);
                ct_index = i < degree/2 ? ct_index : ct_index + degree/2;
                int col_index = (i + j*sq_rt) % lwe_sk_len;
                lwe_ct_tmp[i] = lwe_ct_list[ct_index].a[col_index].ConvertToInt();
            }

            Plaintext lwe_ct_pl;
            batch_encoder.encode(lwe_ct_tmp, lwe_ct_pl);
            evaluator.transform_to_ntt_inplace(lwe_ct_pl, lwe_sk_sqrt_list[j].parms_id());

            if (j == 0) {
                evaluator.multiply_plain(lwe_sk_sqrt_list[j], lwe_ct_pl, result[iter]);
            } else {
                Ciphertext temp;
                evaluator.multiply_plain(lwe_sk_sqrt_list[j], lwe_ct_pl, temp);
                evaluator.add_inplace(result[iter], temp);
            }

        }
    }

    for (int i = 0; i < (int) result.size(); i++) {
        evaluator.transform_from_ntt_inplace(result[i]);
    }

    // sum up all sq_rt tmp results to the first one, each first rotate left one and add to the previous tmp result

    for (int iter = sq_rt-1; iter > 0; iter--) {
        evaluator.rotate_rows_inplace(result[iter], 1, gal_keys);
        evaluator.add_inplace(result[iter-1], result[iter]);
    }

    vector<uint64_t> b_parts(degree);
    for(int i = 0; i < degree; i++){
        b_parts[i] = lwe_ct_list[i].b.ConvertToInt();
    }

    Plaintext lwe_b_pl;
    batch_encoder.encode(b_parts, lwe_b_pl);
    evaluator.negate_inplace(result[0]);
    evaluator.add_plain_inplace(result[0], lwe_b_pl);

    if (gateEval) {
        Plaintext q_shift_pl;
        batch_encoder.encode(q_shift_constant, q_shift_pl);
        evaluator.add_plain_inplace(result[0], q_shift_pl);
    }

    return result[0];
}

vector<vector<int>> generateMatrixU_transpose(int n, const int q = prime_p) {
    vector<vector<int>> U(n,  vector<int>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == 0) {
                U[i][j] = (int) power(primitive_root, j, q);
            } else if (i == n/2) {
                U[i][j] = (int) modInverse(U[0][j], q);
            } else {
                U[i][j] = (int) power(U[i-1][j], 3, q);
            }
        }
    }
    return U;
}

void writeUtemp(const vector<uint64_t> U_temp, const int i) {
    ofstream datafile;
    datafile.open ("../data/U_temp/" + to_string(i) + ".txt");

    for(size_t i = 0; i < U_temp.size(); i++){
        datafile << U_temp[i] << "\n";
    }
    datafile.close();
}

vector<uint64_t> readUtemp(const int i, const int degree=poly_modulus_degree_glb) {
    ifstream datafile;
    datafile.open ("../data/U_temp/" + to_string(i) + ".txt");

    vector<uint64_t> U_temp(degree);
    
    for(size_t i = 0; i < U_temp.size(); i++){
        datafile >> U_temp[i];
    }
    datafile.close();

    return U_temp;
}

// assume that degree/2 has a square root
/**
 * @brief Pre-processed Version.
 * 
 * @param context 
 * @param ct_sqrt_list 
 * @param U_plain_list_c 
 * @param gal_keys 
 * @param degree 
 * @return Ciphertext 
 */
Ciphertext slotToCoeff(const SEALContext& context, const SEALContext& context_coeff, vector<Ciphertext>& ct_sqrt_list, vector<Plaintext>& U_plain_list,
                       const GaloisKeys& gal_keys, const int sq_rt = 128, const int degree=poly_modulus_degree_glb) {
    Evaluator evaluator(context), eval_coeff(context_coeff);
    BatchEncoder batch_encoder(context);

    chrono::high_resolution_clock::time_point time_start, time_end;
    int total_mul = 0;

    vector<Ciphertext> result(sq_rt);
    for (int iter = 0; iter < sq_rt; iter++) {
        for (int j = 0; j < (int) ct_sqrt_list.size(); j++) {
            time_start = chrono::high_resolution_clock::now();
            if (j == 0) {
                evaluator.multiply_plain(ct_sqrt_list[j], U_plain_list[iter * ct_sqrt_list.size() + j], result[iter]);
            } else {
                Ciphertext temp;
                evaluator.multiply_plain(ct_sqrt_list[j], U_plain_list[iter * ct_sqrt_list.size() + j], temp);
                evaluator.add_inplace(result[iter], temp);
            }
            time_end = chrono::high_resolution_clock::now();
            total_mul += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
        }
    }

    for (int i = 0; i < (int) result.size(); i++) {
        evaluator.transform_from_ntt_inplace(result[i]);
    }

    for (int iter = sq_rt-1; iter > 0; iter--) {
        eval_coeff.rotate_rows_inplace(result[iter], 1, gal_keys);
        evaluator.add_inplace(result[iter-1], result[iter]);
    }

    return result[0];
}

inline
long div_mod(long a, long b, long mod = 65537){
    return (a*modInverse(b, mod)) % mod;
}

inline
void mult_scalar_vec(vector<int>& output, const vector<int>& input, int k){
    output.resize(input.size());
    for(size_t i = 0; i < output.size(); i++){
        long temp = ((long)input[i]*(long)k)%65537;
        output[i] = temp;
        if(output[i] < 0)
            cerr <<temp << " " << k << " " << input[i] << endl;
    } 
}

inline
void subtract_two_vec_inplace(vector<int>& output, const vector<int>& input, int numToSolve = -1){
    if(output.size() != input.size())
    {
        cerr << "substracting size not equal." << endl;
    }
    if(numToSolve == -1) numToSolve = input.size();
    for(int i = 0; i < numToSolve; i++){
        output[i] -= input[i];
        output[i] %= 65537; // modulus
        while(output[i] < 0){
            output[i] += 65537;
        }
    }
}

inline
void get_ratio_mult_and_subtract(vector<int>& outputLhs, const vector<int>& inputLhs,
                                 vector<int>& outputRhs, const vector<int>& inputRhs,
                                 const int whichItem, const int numToSolve = -1)
{
    vector<int> temp(inputLhs.size());
    int k = div_mod(outputLhs[whichItem], inputLhs[whichItem]);
    mult_scalar_vec(temp, inputLhs, k);
    subtract_two_vec_inplace(outputLhs, temp);

    mult_scalar_vec(temp, inputRhs, k);
    subtract_two_vec_inplace(outputRhs, temp, numToSolve);
}

inline
vector<long> singleSolve(const long& a, const vector<int>& toSolve, long mod = 65537){
    long a_rev = modInverse(a, mod);
    vector<long> res(toSolve.size());
    for(size_t i = 0; i < toSolve.size(); i++){
        res[i] = ((long)toSolve[i] * a_rev) % 65537;
    }
    return res;
}


vector<vector<long>> equationSolving(vector<vector<int>>& lhs, vector<vector<int>>& rhs, const int& numToSolve = 306){
    vector<int> recoder(lhs[0].size(), -1);
    vector<vector<long>> res(recoder.size());
    size_t counter = 0;

    while(counter < recoder.size()){
        for(size_t i = 0; i < lhs.size(); i++){
            if (lhs[i][counter] != 0){
                if(find(recoder.begin(), recoder.end(), int(i)) != recoder.end()){
                    continue;
                }
                recoder[counter] = i;
                break;
            }
        }

        if(recoder[counter] == -1) {
            // cout << "no solution" << endl;
            return vector<vector<long>>(0);
        }

        for(size_t i = 0; i < lhs.size(); i++) {
	        if ((lhs[i][counter] != 0) && ((int) i != recoder[counter])) {
                get_ratio_mult_and_subtract(lhs[i], lhs[recoder[counter]], rhs[i], rhs[recoder[counter]], counter, numToSolve);
                if (all_of(lhs[i].begin(), lhs[i].end(), [](int j) { return j==0; })) {
                    lhs.erase(lhs.begin() + i);
                    rhs.erase(rhs.begin() + i);
                    return equationSolving(lhs, rhs, numToSolve);
                }
            }
        }
        counter++;
    }

    counter = 0;
    for(size_t i = 0; i < recoder.size(); i++){
        res[i] = singleSolve(lhs[recoder[counter]][counter], rhs[recoder[counter]]);
        counter++;
    }
    return res;
}

// Pick random values to satisfy multi-variable equation.
// For example, given x + y = 10, we might output {2, 8}.
void assignVariable(RandomToStandardAdapter& engine, vector<vector<long>>& res, vector<int>& lhs, int rhs) {
    uniform_int_distribution<uint64_t> dist(0, 65536);

    if (res.size() != lhs.size())
        cerr << "Coefficient and variable size not match." << endl;

    int lastIndex = lhs.size() - 1;

    for (int i = lhs.size(); i > 0; i--) {
        if (lhs[i-1] != 0) {
            lastIndex = i-1;
            break;
        }
    }

    for (int i = 0; i < (int) lhs.size(); i++) {
        if (lhs[i] != 0 && i != lastIndex) {
            res[i][0] = dist(engine);
            long temp = (rhs - (lhs[i] * res[i][0])) % 65537;
            temp = temp < 0 ? temp + 65537 : temp;
            rhs = temp;
        }
    }

    res[lastIndex][0] = div_mod(rhs % 65537, lhs[lastIndex]);
    if (res[lastIndex][0] < 0)
        res[lastIndex][0] += 65537;
}

// Given solved variables with their values, update the remaining equations.
// For example, with equation; x + y + 2z = 10, and z = 2, updated equation would be x + y = 6.
void updateEquation(vector<vector<long>>& res, vector<vector<int>>& lhs, vector<vector<int>>& rhs) {
  for (int i = 0; i < (int) lhs.size(); i++) {
    for (int j = 0; j < (int) res.size(); j++) {
            if (res[j][0] > 0 && lhs[i][j] != 0) {
                long temp = (rhs[i][0] - lhs[i][j] * res[j][0]) % 65537;
                temp = temp < 0 ? temp + 65537 : temp;
                rhs[i][0] = temp;
                lhs[i][j] = 0;
            }
        }
    }
}


vector<vector<long>> equationSolvingRandomBatch(vector<vector<int>>& lhs, vector<vector<int>>& rhs, const int& numToSolve = -1) {
    vector<vector<long>> batched_res(rhs[0].size(), vector<long>(lhs[0].size()));
    prng_seed_type seed;
    for (auto &i : seed) {
        i = random_uint64();
    }

    auto rng = make_shared<Blake2xbPRNGFactory>(Blake2xbPRNGFactory(seed));
    RandomToStandardAdapter engine(rng->create());

    vector<vector<long>> tryRes = equationSolving(lhs, rhs, -1);
    if (tryRes.empty()) {
        for (int k = 0; k < (int) rhs[0].size(); k++) { // for each batched rhs, separate the equation system and solve the variables
            vector<vector<int>> single_rhs(rhs.size(), vector<int>(1)), lhs_copy(lhs.size(), vector<int>(lhs[0].size()));
            for (int i = 0; i < (int) rhs.size(); i++) {
                single_rhs[i][0] = rhs[i][k];
            }
            for (int i = 0; i < (int) lhs.size(); i++) {
                for (int j = 0; j < (int) lhs[0].size(); j++) {
                    lhs_copy[i][j] = lhs[i][j];
                }
            }
            vector<vector<long>> single_res(lhs_copy[0].size(), vector<long>(1));
            while (!lhs_copy.empty()) {
                assignVariable(engine, single_res, lhs_copy[lhs_copy.size() - 1], single_rhs[single_rhs.size() - 1][0]);
                lhs_copy.pop_back();
                single_rhs.pop_back();
                updateEquation(single_res, lhs_copy, single_rhs);
            }

            for (int i = 0; i < (int)lhs[0].size(); i++) {
                batched_res[k][i] = single_res[i][0];
            }
        }
    } else {
      return tryRes;
    }
    return batched_res;
}


vector<vector<int>> generateIdentityMatrix(uint64_t size) {
  vector<vector<int>> identity(size);

  for (int i = 0; i < (int) size; i++) {
    identity[i].resize(size);
    for (int j = 0; j < (int) size; j++) {
      identity[i][j] = (i == j);
    }
  }

  return identity;
}

vector<vector<uint64_t>> findInverse(vector<vector<int>> mat, uint64_t q) {
  vector<vector<int>> identity = generateIdentityMatrix(mat.size());
  
  vector<vector<long>> result = equationSolvingRandomBatch(mat, identity, -1);

  vector<vector<uint64_t>> ret((int)result.size());
  for (int i = 0; i < (int) result.size(); i++) {
    ret[i].resize(result[i].size());
    for (int j = 0; j < (int) result[0].size(); j++) {
      ret[i][j] = (uint64_t) result[i][j];
    }
  }

  return ret;
  
}

Ciphertext inverse_slotToCoeff_WOPrepreocess(const SEALContext& context, const SEALContext& context_coeff,
					     vector<Ciphertext>& ct_sqrt_list, const GaloisKeys& gal_keys,
					     const int degree=poly_modulus_degree_glb, const int q=65537,
					     const uint64_t scalar=1) {
    Evaluator evaluator(context), eval_coeff(context_coeff);
    BatchEncoder batch_encoder(context);

    chrono::high_resolution_clock::time_point time_start, time_end;
    uint64_t total_U = 0;

    time_start = chrono::high_resolution_clock::now();
    vector<vector<int>> U_old = generateMatrixU_transpose(degree, q);
    vector<vector<uint64_t>> U_inverse = findInverse(U_old, q);
    time_end = chrono::high_resolution_clock::now();
    total_U += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    int sq_rt = sqrt(degree/2);


    vector<Ciphertext> result(sq_rt);
    for (int iter = 0; iter < sq_rt; iter++) {
        for (int j = 0; j < (int) ct_sqrt_list.size(); j++) {

            time_start = chrono::high_resolution_clock::now();
            vector<uint64_t> U_tmp(degree);
            for (int i = 0; i < degree; i++) {
                int row_index = (i-iter) % (degree/2) < 0 ? (i-iter) % (degree/2) + degree/2 : (i-iter) % (degree/2);
                row_index = i < degree/2 ? row_index : row_index + degree/2;
                int col_index = (i + j*sq_rt) % (degree/2);
                if (j < (int) ct_sqrt_list.size() / 2) { // first half
                    col_index = i < degree/2 ? col_index : col_index + degree/2;
                } else {
                    col_index = i < degree/2 ? col_index + degree/2 : col_index;
                }
                U_tmp[i] = (U_inverse[row_index][col_index] * scalar) % q;
            }

            Plaintext U_plain;
            batch_encoder.encode(U_tmp, U_plain);
            evaluator.transform_to_ntt_inplace(U_plain, ct_sqrt_list[j].parms_id());

            time_end = chrono::high_resolution_clock::now();
            total_U += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

            if (j == 0) {
                evaluator.multiply_plain(ct_sqrt_list[j], U_plain, result[iter]);
            } else {
                Ciphertext temp;
                evaluator.multiply_plain(ct_sqrt_list[j], U_plain, temp);
                evaluator.add_inplace(result[iter], temp);
            }
        }
    }

    for (int i = 0; i < (int) result.size(); i++) {
        evaluator.transform_from_ntt_inplace(result[i]);
    }

    for (int iter = sq_rt-1; iter > 0; iter--) {
        eval_coeff.rotate_rows_inplace(result[iter], 1, gal_keys);
        evaluator.add_inplace(result[iter-1], result[iter]);
    }

    cout << "    TOTAL process U time: " << total_U << endl;

    return result[0];
}

/**
 * @brief Non Pre-processed Version.
 * 
 * @param context 
 * @param ct_sqrt_list 
 * @param U_plain_list_c 
 * @param gal_keys 
 * @param degree 
 * @return Ciphertext 
 */
Ciphertext slotToCoeff_WOPrepreocess(const SEALContext& context, const SEALContext& context_coeff,
				     vector<Ciphertext>& ct_sqrt_list, const GaloisKeys& gal_keys,
                                     const int degree=poly_modulus_degree_glb, const int q = prime_p,
				     const uint64_t scalar=1) {
    Evaluator evaluator(context), eval_coeff(context_coeff);
    BatchEncoder batch_encoder(context);

    chrono::high_resolution_clock::time_point time_start, time_end;
    uint64_t total_U = 0;

    time_start = chrono::high_resolution_clock::now();
    vector<vector<int>> U = generateMatrixU_transpose(degree, q);
    time_end = chrono::high_resolution_clock::now();
    total_U += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    int sq_rt = sqrt(degree/2);

    vector<Ciphertext> result(sq_rt);
    for (int iter = 0; iter < sq_rt; iter++) {
        for (int j = 0; j < (int) ct_sqrt_list.size(); j++) {

            time_start = chrono::high_resolution_clock::now();
            vector<uint64_t> U_tmp(degree);
            for (int i = 0; i < degree; i++) {
                int row_index = (i-iter) % (degree/2) < 0 ? (i-iter) % (degree/2) + degree/2 : (i-iter) % (degree/2);
                row_index = i < degree/2 ? row_index : row_index + degree/2;
                int col_index = (i + j*sq_rt) % (degree/2);
                if (j < (int) ct_sqrt_list.size() / 2) { // first half
                    col_index = i < degree/2 ? col_index : col_index + degree/2;
                } else {
                    col_index = i < degree/2 ? col_index + degree/2 : col_index;
                }
                U_tmp[i] = (U[row_index][col_index] * scalar) % q;
            }
            writeUtemp(U_tmp, j*sq_rt + iter);

            Plaintext U_plain;
            batch_encoder.encode(U_tmp, U_plain);
            evaluator.transform_to_ntt_inplace(U_plain, ct_sqrt_list[j].parms_id());

            time_end = chrono::high_resolution_clock::now();
            total_U += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

            if (j == 0) {
                evaluator.multiply_plain(ct_sqrt_list[j], U_plain, result[iter]);
            } else {
                Ciphertext temp;
                evaluator.multiply_plain(ct_sqrt_list[j], U_plain, temp);
                evaluator.add_inplace(result[iter], temp);
            }
        }
    }

    for (int i = 0; i < (int) result.size(); i++) {
        evaluator.transform_from_ntt_inplace(result[i]);
    }

    for (int iter = sq_rt-1; iter > 0; iter--) {
        eval_coeff.rotate_rows_inplace(result[iter], 1, gal_keys);
        evaluator.add_inplace(result[iter-1], result[iter]);
    }

    cout << "    TOTAL process U time: " << total_U << endl;

    return result[0];
}

// bootstrap rangecheck function that calculate a poly given error bound, and no care for input outside domain
void Bootstrap_FastRangeCheck_Random(const SecretKey& bfv_secret_key, Ciphertext& output, const Ciphertext& input, const size_t& degree, const RelinKeys &relin_keys,
                                     const SEALContext& context, const vector<uint64_t>& rangeCheckIndices, const int firstLevel, const int secondLevel,
                                     map<int, bool>& firstLevelMod, map<int, bool>& secondLevelMod, const int f_zero = 0) {
    MemoryPoolHandle my_pool_larger = MemoryPoolHandle::New(true);
    auto old_prof_larger = MemoryManager::SwitchProfile(std::make_unique<MMProfFixed>(std::move(my_pool_larger)));

    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);
    Decryptor decryptor(context, bfv_secret_key);
    vector<Ciphertext> kCTs(firstLevel), kToMCTs(secondLevel);

    calUptoDegreeK_bigPrime(kCTs, input, firstLevel, relin_keys, context, firstLevelMod);
    cout << "   Noise after first level: " << decryptor.invariant_noise_budget(kCTs[kCTs.size()-1]) << " bits\n";
    calUptoDegreeK_bigPrime(kToMCTs, kCTs[kCTs.size()-1], secondLevel, relin_keys, context, secondLevelMod);
    cout << "   Noise after second level: " << decryptor.invariant_noise_budget(kToMCTs[kToMCTs.size()-1]) << " bits\n";

    // for (int j = 0; j < (int) kCTs.size(); j++) {
    //     evaluator.mod_switch_to_inplace(kCTs[j], kToMCTs[kToMCTs.size()-1].parms_id());
    // }
    // for (int j = 0; j < (int) kToMCTs.size(); j++) {
    //     evaluator.mod_switch_to_next_inplace(kToMCTs[j]);
    // }

    Ciphertext temp_relin(context);

    Plaintext plainInd;
    plainInd.resize(degree);
    plainInd.parms_id() = parms_id_zero;
    for (int i = 0; i < (int) degree; i++) {
        plainInd.data()[i] = 0;
    }

    for(int i = 0; i < secondLevel; i++) {
        Ciphertext levelSum;
        bool flag = false;
        for(int j = 0; j < firstLevel; j++) {
            if(rangeCheckIndices[i*firstLevel+j] != 0) {
                plainInd.data()[0] = rangeCheckIndices[i*firstLevel+j];
                if (!flag) {
                    evaluator.multiply_plain(kCTs[j], plainInd, levelSum);
                    flag = true;
                } else {
                    Ciphertext tmp;
                    evaluator.multiply_plain(kCTs[j], plainInd, tmp);
                    evaluator.add_inplace(levelSum, tmp);
                }
            }
        }
        evaluator.mod_switch_to_inplace(levelSum, kToMCTs[i].parms_id()); // mod down the plaintext multiplication noise
        if(i == 0) {
            output = levelSum;
        } else if (i == 1) { // initialize for temp_relin, which is of ct size = 3
            evaluator.multiply(levelSum, kToMCTs[i - 1], temp_relin);
        } else {
            evaluator.multiply_inplace(levelSum, kToMCTs[i - 1]);
            evaluator.add_inplace(temp_relin, levelSum);
        }
    }

    for(int i = 0; i < firstLevel; i++){
        kCTs[i].release();
    }
    for(int i = 0; i < secondLevel; i++){
        kToMCTs[i].release();
    }

    evaluator.relinearize_inplace(temp_relin, relin_keys);
    evaluator.add_inplace(output, temp_relin);
    temp_relin.release();

    if (f_zero) {
        plainInd.data()[0] = f_zero;
        evaluator.add_plain_inplace(output, plainInd);
    }

    MemoryManager::SwitchProfile(std::move(old_prof_larger));
}

Ciphertext calculateDegree(const SEALContext& context, const RelinKeys &relin_keys, Ciphertext& input, map<int, bool> modDownIndices, int degree) {
    Evaluator evaluator(context);

    vector<Ciphertext> output(degree);
    output[0] = input;
    vector<int> calculated(degree, 0), numMod(degree, 0);
    calculated[0] = 1;

    Ciphertext res, base;

    auto toCalculate = degree;
    int resdeg = 0;
    int basedeg = 1;
    base = input;
    while(toCalculate > 0){
        if(toCalculate & 1){
            toCalculate -= 1;
            resdeg += basedeg;
            if(calculated[resdeg-1] != 0){
                res = output[resdeg - 1];
            } else {
                if(resdeg == basedeg){
                    res = base; // should've never be used as base should have made calculated[basedeg-1]
                } else {
                    numMod[resdeg-1] = numMod[basedeg-1];

                    evaluator.mod_switch_to_inplace(res, base.parms_id()); // match modulus
                    evaluator.multiply_inplace(res, base);
                    evaluator.relinearize_inplace(res, relin_keys);
                    if(modDownIndices.count(resdeg) && !modDownIndices[resdeg]) {
                        modDownIndices[resdeg] = true;
                        evaluator.mod_switch_to_next_inplace(res);
                        numMod[resdeg-1]+=1;
                    }
                }
                output[resdeg-1] = res;
                calculated[resdeg-1] += 1;
            }
        } else {
            toCalculate /= 2;
            basedeg *= 2;
            if(calculated[basedeg-1] != 0){
                base = output[basedeg - 1];
            } else {
                numMod[basedeg-1] = numMod[basedeg/2-1];
                evaluator.square_inplace(base);
                evaluator.relinearize_inplace(base, relin_keys);
                while(modDownIndices.count(basedeg) && !modDownIndices[basedeg]) {
                    modDownIndices[basedeg] = true;
                    evaluator.mod_switch_to_next_inplace(base);
                    numMod[basedeg-1]+=1;
                }
                output[basedeg-1] = base;
                calculated[basedeg-1] += 1;
            }
        }
    }

    return output[output.size()-1];
}


Ciphertext raisePowerToPrime(const SEALContext& context, const RelinKeys &relin_keys, Ciphertext& input, map<int, bool> modDownIndices_1, map<int, bool> modDownIndices_2,
                             int degree_1, int degree_2, int prime = prime_p) {

    Ciphertext tmp = calculateDegree(context, relin_keys, input, modDownIndices_1, degree_1);
    tmp = calculateDegree(context, relin_keys, tmp, modDownIndices_2, degree_2);

    return tmp;
}


// bootstrap rangecheck function that calculate a poly given error bound, condition mapping, raise random result to 1
void Bootstrap_FastRangeCheck_Condition(SecretKey& bfv_secret_key, Ciphertext& output, const Ciphertext& input, const size_t& degree, const RelinKeys &relin_keys,
                                        const SEALContext& context, const vector<uint64_t>& rangeCheckIndices, const int eval_level1, const int eval_level2,
                                        map<int, bool>& eval_mod1, map<int, bool>& eval_mod2, const int raise_level1, const int raise_level2,
                                        map<int, bool>& raise_mod1, map<int, bool>& raise_mod2) {
    MemoryPoolHandle my_pool_larger = MemoryPoolHandle::New(true);
    auto old_prof_larger = MemoryManager::SwitchProfile(std::make_unique<MMProfFixed>(std::move(my_pool_larger)));

    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);
    Decryptor decryptor(context, bfv_secret_key);

    chrono::high_resolution_clock::time_point time_start, time_end;
    time_start = chrono::high_resolution_clock::now();

    vector<Ciphertext> kCTs(eval_level1), kToMCTs(eval_level2);

    cout << "first: " << eval_level1 << ", second: " << eval_level2 << endl; 

    calUptoDegreeK_bigPrime(kCTs, input, eval_level1, relin_keys, context, eval_mod1);
    calUptoDegreeK_bigPrime(kToMCTs, kCTs[kCTs.size()-1], eval_level2, relin_keys, context, eval_mod2);

    for (int j = 0; j < (int) kCTs.size(); j++) {
        evaluator.mod_switch_to_inplace(kCTs[j], kToMCTs[kToMCTs.size()-1].parms_id());
    }
    // for (int j = 0; j < (int) kToMCTs.size(); j++) {
    //     evaluator.mod_switch_to_next_inplace(kToMCTs[j]);
    // }

    Ciphertext temp_relin(context);

    Plaintext plainInd;
    plainInd.resize(degree);
    plainInd.parms_id() = parms_id_zero;
    for (int i = 0; i < (int) degree; i++) {
        plainInd.data()[i] = 0;
    }

    for(int i = 0; i < eval_level2; i++) {
        Ciphertext levelSum;
        bool flag = false;
        for(int j = 0; j < eval_level1; j++) {
            if(rangeCheckIndices[i*eval_level1+j] != 0) {
                plainInd.data()[0] = rangeCheckIndices[i*eval_level1+j];
                if (!flag) {
                    evaluator.multiply_plain(kCTs[j], plainInd, levelSum);
                    flag = true;
                } else {
                    Ciphertext tmp;
                    evaluator.multiply_plain(kCTs[j], plainInd, tmp);
                    evaluator.add_inplace(levelSum, tmp);
                }
            }
        }
        evaluator.mod_switch_to_inplace(levelSum, kToMCTs[i].parms_id()); // mod down the plaintext multiplication noise
        if(i == 0) {
            output = levelSum;
        } else if (i == 1) { // initialize for temp_relin, which is of ct size = 3
            evaluator.multiply(levelSum, kToMCTs[i - 1], temp_relin);
        } else {
            evaluator.multiply_inplace(levelSum, kToMCTs[i - 1]);
            evaluator.add_inplace(temp_relin, levelSum);
        }
    }

    for(int i = 0; i < eval_level1; i++){
        kCTs[i].release();
    }
    for(int i = 0; i < eval_level2; i++){
        kToMCTs[i].release();
    }

    evaluator.relinearize_inplace(temp_relin, relin_keys);
    evaluator.add_inplace(output, temp_relin);
    temp_relin.release();

    time_end = chrono::high_resolution_clock::now();
    cout << "   first evaluation half: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;

    time_start = chrono::high_resolution_clock::now();
    output = raisePowerToPrime(context, relin_keys, output, raise_mod1, raise_mod2, raise_level1, raise_level2, prime_p);
    time_end = chrono::high_resolution_clock::now();
    cout << "   second raise power half: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;

    MemoryManager::SwitchProfile(std::move(old_prof_larger));
}



void Bootstrap_RangeCheck_PatersonStockmeyer_bigPrime(Ciphertext& ciphertext, const Ciphertext& input, const vector<uint64_t>& rangeCheckIndices,
                                             const int modulus, const size_t& degree, const RelinKeys &relin_keys, const SEALContext& context, const SecretKey& bfv_secret_key,
                                             const int f_zero = 0, const bool gateEval = false, const bool skip_first_odd = false,
                                             const int firstDegree = 256, const int secondDegree = 256) {
    MemoryPoolHandle my_pool_larger = MemoryPoolHandle::New(true);
    auto old_prof_larger = MemoryManager::SwitchProfile(std::make_unique<MMProfFixed>(std::move(my_pool_larger)));

    Decryptor decryptor(context, bfv_secret_key);
    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);
    vector<Ciphertext> kCTs(firstDegree);

    map<int, bool> modDownIndices_firstLevel = {{2, false}, {8, false}, {32, false}, {128, false}, {512, false}};
    map<int, bool> modDownIndices_secondLevel = {{2, false}, {8, false}, {32, false}, {64, false}, {128, false}, {512, false}};

    calUptoDegreeK_bigPrime(kCTs, input, firstDegree, relin_keys, context, modDownIndices_firstLevel, skip_first_odd);

    vector<Ciphertext> kToMCTs(secondDegree);
    calUptoDegreeK_bigPrime(kToMCTs, kCTs[kCTs.size()-1], secondDegree, relin_keys, context, modDownIndices_secondLevel);

    // for (int j = 0; j < (int) kToMCTs.size(); j++) {
    //     evaluator.mod_switch_to_next_inplace(kToMCTs[j]);
    // }
    for (int j = 0; j < (int) kCTs.size(); j++) {
        evaluator.mod_switch_to_inplace(kCTs[j], kToMCTs[kToMCTs.size()-1].parms_id());
    }
    // cout << "Noise for last: " << decryptor.invariant_noise_budget(kToMCTs[kToMCTs.size()-1]) << " bits\n";
    for (int j = 0; j < (int) kToMCTs.size(); j++) {
        evaluator.mod_switch_to_next_inplace(kToMCTs[j]);
    }
    // cout << "Noise for last: " << decryptor.invariant_noise_budget(kToMCTs[kToMCTs.size()-1]) << " bits\n";
    

    Ciphertext temp_relin(context);

    Plaintext plainInd;
    plainInd.resize(degree);
    plainInd.parms_id() = parms_id_zero;
    for (int i = 0; i < (int) degree; i++) {
        plainInd.data()[i] = 0;
    }

    for(int i = 0; i < secondDegree; i++) {
        Ciphertext levelSum;
        bool flag = false;
        for(int j = 0; j < firstDegree; j++) {
            if(rangeCheckIndices[i*firstDegree+j] != 0) {
                plainInd.data()[0] = rangeCheckIndices[i*firstDegree+j];
                if (!flag) {
                    evaluator.multiply_plain(kCTs[j], plainInd, levelSum);
                    flag = true;
                } else {
                    Ciphertext tmp;
                    evaluator.multiply_plain(kCTs[j], plainInd, tmp);
                    evaluator.add_inplace(levelSum, tmp);
                }
            }
        }
        evaluator.mod_switch_to_inplace(levelSum, kToMCTs[i].parms_id()); // mod down the plaintext multiplication noise
        if(i == 0) {
            ciphertext = levelSum;
        } else if (i == 1) { // initialize for temp_relin, which is of ct size = 3
            evaluator.multiply(levelSum, kToMCTs[i - 1], temp_relin);
        } else {
            evaluator.multiply_inplace(levelSum, kToMCTs[i - 1]);
            evaluator.add_inplace(temp_relin, levelSum);
        }
    }

    for(int i = 0; i < firstDegree; i++){
        kCTs[i].release();
    }
    for(int i = 0; i < secondDegree; i++){
        kToMCTs[i].release();
    }

    evaluator.relinearize_inplace(temp_relin, relin_keys);
    evaluator.add_inplace(ciphertext, temp_relin);
    temp_relin.release();

    plainInd.data()[0] = f_zero;
    evaluator.negate_inplace(ciphertext);
    evaluator.add_plain_inplace(ciphertext, plainInd);

    cout << "Noise after function: " << decryptor.invariant_noise_budget(ciphertext) << " bits\n";

    MemoryManager::SwitchProfile(std::move(old_prof_larger));
}


void Bootstrap_RangeCheck_PatersonStockmeyer(Ciphertext& ciphertext, const Ciphertext& input, const vector<uint64_t>& rangeCheckIndices,
                                             const int modulus, const size_t& degree, const RelinKeys &relin_keys, const SEALContext& context, const SecretKey& bfv_secret_key, 
                                             const int f_zero = 0, const bool gateEval = false, const bool skip_first_odd = true,
                                             const int firstDegree = 256, const int secondDegree = 256) {
    MemoryPoolHandle my_pool_larger = MemoryPoolHandle::New(true);
    auto old_prof_larger = MemoryManager::SwitchProfile(std::make_unique<MMProfFixed>(std::move(my_pool_larger)));

    Decryptor decryptor(context, bfv_secret_key);
    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);
    vector<Ciphertext> kCTs(firstDegree), kToMCTs(secondDegree);

    calUptoDegreeK(kCTs, input, firstDegree, relin_keys, context, skip_first_odd);
    calUptoDegreeK(kToMCTs, kCTs[kCTs.size()-1], secondDegree, relin_keys, context);

    for (int j = 0; j < (int) kCTs.size(); j++) {
        evaluator.mod_switch_to_inplace(kCTs[j], kToMCTs[kToMCTs.size()-1].parms_id());
    }
    for (int j = 0; j < (int) kToMCTs.size(); j++) {
        evaluator.mod_switch_to_next_inplace(kToMCTs[j]);
    }

    Ciphertext temp_relin(context);
    
    Plaintext plainInd;
    plainInd.resize(degree);
    plainInd.parms_id() = parms_id_zero;
    for (int i = 0; i < (int) degree; i++) {
        plainInd.data()[i] = 0;
    }

    for(int i = 0; i < secondDegree; i++) {
        Ciphertext levelSum;
        bool flag = false;
        for(int j = 0; j < firstDegree; j++) {
            if(rangeCheckIndices[i*firstDegree+j] != 0) {
                plainInd.data()[0] = rangeCheckIndices[i*firstDegree+j];
                if (!flag) {
                    evaluator.multiply_plain(kCTs[j], plainInd, levelSum);
                    flag = true;
                } else {
                    Ciphertext tmp;
                    evaluator.multiply_plain(kCTs[j], plainInd, tmp);
                    evaluator.add_inplace(levelSum, tmp);
                }
            }
        }
        evaluator.mod_switch_to_inplace(levelSum, kToMCTs[i].parms_id()); // mod down the plaintext multiplication noise
        if(i == 0) {
            ciphertext = levelSum;
        } else if (i == 1) { // initialize for temp_relin, which is of ct size = 3
            evaluator.multiply(levelSum, kToMCTs[i - 1], temp_relin);
        } else {
            evaluator.multiply_inplace(levelSum, kToMCTs[i - 1]);
            evaluator.add_inplace(temp_relin, levelSum);
        }
    }
    
    for(int i = 0; i < firstDegree; i++){
        kCTs[i].release();
    }
    for(int i = 0; i < secondDegree; i++){
        kToMCTs[i].release();
    }

    evaluator.relinearize_inplace(temp_relin, relin_keys);
    evaluator.add_inplace(ciphertext, temp_relin);
    temp_relin.release();


    plainInd.data()[0] = f_zero;
    evaluator.negate_inplace(ciphertext);
    evaluator.add_plain_inplace(ciphertext, plainInd);

    cout << "Noise after function: " << decryptor.invariant_noise_budget(ciphertext) << " bits\n";

    if (gateEval) { // flip 0 to q/3, q/3 to 0
        plainInd.data()[0] = modulus/3;
        evaluator.negate_inplace(ciphertext);
        evaluator.add_plain_inplace(ciphertext, plainInd);
    }
    
    MemoryManager::SwitchProfile(std::move(old_prof_larger));
}


Ciphertext encryptLWEskUnderBFV(const SEALContext& context, const size_t& degree,
                                const PublicKey& BFVpk, const SecretKey& BFVsk,
                                const regevSK& regSk, const regevParam& params) { 
    Ciphertext switchingKey(context);

    BatchEncoder batch_encoder(context);
    Encryptor encryptor(context, BFVpk);
    encryptor.set_secret_key(BFVsk);

    int tempn = 1;
    for(tempn = 1; tempn < params.n; tempn *= 2){}

    vector<uint64_t> skInt(degree);
    for(size_t i = 0; i < degree; i++){
        auto tempindex = i%uint64_t(tempn);
        if(int(tempindex) >= params.n) {
            skInt[i] = 0;
        } else {
            skInt[i] = uint64_t(regSk[tempindex].ConvertToInt() % params.q);
        }
    }
    Plaintext plaintext;
    batch_encoder.encode(skInt, plaintext);
    encryptor.encrypt_symmetric(plaintext, switchingKey);

    return switchingKey;
}


vector<regevCiphertext> extractRLWECiphertextToLWECiphertext(Ciphertext& rlwe_ct, const int ring_dim = poly_modulus_degree_glb,
                                                             const int n = 1024, const int p = prime_p, const uint64_t big_prime = 1152921504589938689) {
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


vector<regevCiphertext> preprocess_NAND(const vector<regevCiphertext>& ct_list_1, const vector<regevCiphertext>& ct_list_2,
                                        const regevParam& params) {

    chrono::high_resolution_clock::time_point time_start, time_end;
    time_start = chrono::high_resolution_clock::now();

    vector<regevCiphertext> result(ct_list_1.size());
    for (int i = 0; i < (int) ct_list_1.size(); i++) {
        result[i].a = NativeVector(params.n);
        for (int j = 0; j < params.n; j++) {
            result[i].a[j] = (ct_list_1[i].a[j].ConvertToInt() + ct_list_2[i].a[j].ConvertToInt()) % params.q;
        }
        result[i].b = (ct_list_1[i].b.ConvertToInt() + ct_list_2[i].b.ConvertToInt()) % params.q;
    }

    time_end = chrono::high_resolution_clock::now();
    cout << "TOTAL prepare NAND input time: " << chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    
    return result;
}


void modDownToPrime(Ciphertext& coeff, const int ring_dim, const uint64_t big_prime, const uint64_t small_prime=268369921) {
    for (int i = 0; i < ring_dim; i++) {
        coeff.data(1)[i] =  uint64_t(float(coeff.data(1)[i]) * float(small_prime) / float(big_prime));
        coeff.data(0)[i] =  uint64_t(float(coeff.data(0)[i]) * float(small_prime) / float(big_prime));
    }
}


vector<regevCiphertext> bootstrap_bigPrime(vector<regevCiphertext>& lwe_ct_list, Ciphertext& lwe_sk_encrypted, const SEALContext& seal_context,
                                           const SEALContext& seal_context_last, const RelinKeys& relin_keys, const GaloisKeys& gal_keys,
                                           const GaloisKeys& gal_keys_coeff, const int ring_dim, const int n, const int p, const KSwitchKeys& ksk,
                                           const vector<uint64_t>& rangeCheckIndices, const MemoryPoolHandle& my_pool, const SecretKey& bfv_secret_key,
                                           const vector<uint64_t>& q_shift_constant, const int f_zero = 0, const bool gateEval = false,
                                           const bool skip_first_odd = true, const int firstDegree = 256, const int secondDegree = 256,
                                           const uint64_t bigPrime = 1152921504581877761, const uint64_t smallPrime = 268369921) {
    chrono::high_resolution_clock::time_point time_start, time_end;
    uint64_t total_preprocess = 0, total_online = 0;

    Evaluator evaluator(seal_context);
    BatchEncoder batch_encoder(seal_context);
    Decryptor decryptor(seal_context, bfv_secret_key);

    int sq_sk = sqrt(n), sq_ct = sqrt(ring_dim/2);
    vector<Ciphertext> lwe_sk_sqrt_list(sq_sk), ct_sqrt_list(2*sq_ct);

    Ciphertext lwe_sk_column;

    time_start = chrono::high_resolution_clock::now();
    evaluator.rotate_columns(lwe_sk_encrypted, gal_keys, lwe_sk_column);
    for (int i = 0; i < sq_sk; i++) {
        evaluator.rotate_rows(lwe_sk_encrypted, sq_sk * i, gal_keys, lwe_sk_sqrt_list[i]);
        evaluator.transform_to_ntt_inplace(lwe_sk_sqrt_list[i]);
    }
    time_end = chrono::high_resolution_clock::now();
    total_preprocess += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

    time_start = chrono::high_resolution_clock::now();
    Ciphertext result = evaluatePackedLWECiphertext(seal_context, lwe_ct_list, lwe_sk_sqrt_list, gal_keys, n, q_shift_constant, ring_dim, gateEval);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    cout << "TOTAL TIME for evaluation: " << total_online << endl;
    // cout << "Noise: " << decryptor.invariant_noise_budget(result) << " bits\n";

    // evaluator.mod_switch_to_next_inplace(result);
    // cout << "Noise after mod down: " << decryptor.invariant_noise_budget(result) << " bits\n";

    Plaintext pl;
    vector<uint64_t> v(ring_dim);
    // decryptor.decrypt(result, pl);
    // batch_encoder.decode(pl, v);
    // cout << "Decrypt after evaluation: \n" << v << endl;

    Ciphertext range_check_res;
    time_start = chrono::high_resolution_clock::now();
    Bootstrap_RangeCheck_PatersonStockmeyer_bigPrime(range_check_res, result, rangeCheckIndices, p, ring_dim,
                                                     relin_keys, seal_context, bfv_secret_key, f_zero, gateEval, skip_first_odd,
                                                     firstDegree, secondDegree);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    cout << "TOTAL TIME for rangecheck: " << total_online << endl;

    // decryptor.decrypt(range_check_res, pl);
    // batch_encoder.decode(pl, v);
    // cout << "Decrypt after rangeCheck: \n" << v << endl;

    ////////////////////////////////////////// SLOT TO COEFFICIENT /////////////////////////////////////////////////////

    time_start = chrono::high_resolution_clock::now();
    evaluator.mod_switch_to_next_inplace(range_check_res);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    // cout << "Noise after range check: " << decryptor.invariant_noise_budget(range_check_res) << " bits\n";

    time_start = chrono::high_resolution_clock::now();
    Ciphertext range_check_res_copy(range_check_res);

    Evaluator eval_coeff(seal_context_last);
    eval_coeff.rotate_columns_inplace(range_check_res_copy, gal_keys_coeff);
    for (int i = 0; i < sq_ct; i++) {
        eval_coeff.rotate_rows(range_check_res, sq_ct * i, gal_keys_coeff, ct_sqrt_list[i]);
        eval_coeff.transform_to_ntt_inplace(ct_sqrt_list[i]);
        eval_coeff.rotate_rows(range_check_res_copy, sq_ct * i, gal_keys_coeff, ct_sqrt_list[i+sq_ct]);
        eval_coeff.transform_to_ntt_inplace(ct_sqrt_list[i+sq_ct]);
    }

    // vector<Plaintext> U_plain_list(ring_dim);
    // for (int iter = 0; iter < sq_ct; iter++) {
    //     for (int j = 0; j < (int) ct_sqrt_list.size(); j++) {
    //         vector<uint64_t> U_tmp = readUtemp(j*sq_ct + iter);
    //         batch_encoder.encode(U_tmp, U_plain_list[iter * ct_sqrt_list.size() + j]);
    //         evaluator.transform_to_ntt_inplace(U_plain_list[iter * ct_sqrt_list.size() + j], ct_sqrt_list[j].parms_id());
    //     }
    // }
    time_end = chrono::high_resolution_clock::now();
    total_preprocess += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

    time_start = chrono::high_resolution_clock::now();
    // Ciphertext coeff = slotToCoeff(seal_context, seal_context_last, ct_sqrt_list, U_plain_list, gal_keys_coeff, 128, ring_dim);
    Ciphertext coeff = slotToCoeff_WOPrepreocess(seal_context, seal_context_last, ct_sqrt_list, gal_keys_coeff, 128, ring_dim, p);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    cout << "TOTAL TIME for slotToCoeff: " << total_online << endl;

    // cout << "PLAINTEX OF SLOTTOCOEFF\n";
    // decryptor.decrypt(coeff, pl);
    // for (int i = 0; i < ring_dim; i++) {
    //   cout << pl[i] << ",";
    // }
    // cout << endl;


    ////////////////////////////////////////////////// MANUAL MOD DOWN /////////////////////////////////////////////////

    // modDownToPrime(coeff, ring_dim, bigPrime, smallPrime);


    ////////////////////////////////////////////////// KEY SWITCHING ///////////////////////////////////////////////////

    time_start = chrono::high_resolution_clock::now();
    while(seal_context.last_parms_id() != coeff.parms_id()){
        evaluator.mod_switch_to_next_inplace(coeff);
    }
    // cout << "Noise before key switch: " << decryptor.invariant_noise_budget(coeff) << " bits\n";

    Ciphertext copy_coeff = coeff;
    auto ct_in_iter = util::iter(copy_coeff);
    ct_in_iter += coeff.size() - 1;
    seal::util::set_zero_poly(ring_dim, 1, coeff.data(1)); // notice that the coeff_mod.size() is hardcoded to 1, thus this needs to be performed on the last level

    evaluator.switch_key_inplace(coeff, *ct_in_iter, static_cast<const KSwitchKeys &>(ksk), 0, my_pool);

    // cout << "Noise before extraction: " << decryptor.invariant_noise_budget(coeff) << " bits\n";

    vector<regevCiphertext> lwe_ct_results = extractRLWECiphertextToLWECiphertext(coeff, ring_dim, n, p, bigPrime);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    cout << "TOTAL TIME for Extraction: " << total_online << endl;
    cout << "TOTAL PREPROCESS TIME: " << total_preprocess << endl;
    cout << "TOTAL ONLINE TIME: " << total_online << endl;

    return lwe_ct_results;
}

vector<regevCiphertext> bootstrap(vector<regevCiphertext>& lwe_ct_list, Ciphertext& lwe_sk_encrypted, const SEALContext& seal_context,
                                  const SEALContext& seal_context_last, const RelinKeys& relin_keys, const GaloisKeys& gal_keys,
                                  const GaloisKeys& gal_keys_coeff, const int ring_dim, const int n, const int p, const KSwitchKeys& ksk,
                                  const vector<uint64_t>& rangeCheckIndices, const MemoryPoolHandle& my_pool, const SecretKey& bfv_secret_key,
                                  const vector<uint64_t>& q_shift_constant, const int f_zero = 0, const bool gateEval = false,
                                  const bool skip_first_odd = true, const int firstDegree = 256, const int secondDegree = 256,
                                  const int sq_ct = 128, const int sq_rt = 128) {
    chrono::high_resolution_clock::time_point time_start, time_end;
    uint64_t total_preprocess = 0, total_online = 0;

    Evaluator evaluator(seal_context);
    BatchEncoder batch_encoder(seal_context);
    Decryptor decryptor(seal_context, bfv_secret_key);

    int sq_sk = sqrt(n);
    vector<Ciphertext> lwe_sk_sqrt_list(sq_sk), ct_sqrt_list(2*sq_ct);

    Ciphertext lwe_sk_column;

    time_start = chrono::high_resolution_clock::now();

    evaluator.rotate_columns(lwe_sk_encrypted, gal_keys, lwe_sk_column);
    for (int i = 0; i < sq_sk; i++) {
        evaluator.rotate_rows(lwe_sk_encrypted, sq_sk * i, gal_keys, lwe_sk_sqrt_list[i]);
        evaluator.transform_to_ntt_inplace(lwe_sk_sqrt_list[i]);
    }
    time_end = chrono::high_resolution_clock::now();
    total_preprocess += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

    time_start = chrono::high_resolution_clock::now();
    Ciphertext result = evaluatePackedLWECiphertext(seal_context, lwe_ct_list, lwe_sk_sqrt_list, gal_keys, n, q_shift_constant, ring_dim, gateEval);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    cout << "TOTAL TIME for evaluation: " << total_online << endl;
    cout << "Noise: " << decryptor.invariant_noise_budget(result) << " bits\n";

    Plaintext pl;
    vector<uint64_t> v(ring_dim);
    // decryptor.decrypt(result, pl);
    // batch_encoder.decode(pl, v);
    // cout << "Decrypt after evaluation: \n" << v << endl;

    Ciphertext range_check_res;
    time_start = chrono::high_resolution_clock::now();
    if (gateEval) {
        map<int, bool> modDownIndices1 = {{4, false}};
        map<int, bool> modDownIndices2 = {{4, false}, {16, false}};
        Bootstrap_FastRangeCheck_Random(bfv_secret_key, range_check_res, result, ring_dim, relin_keys, seal_context, rangeCheckIndices,
                                        firstDegree, secondDegree, modDownIndices1, modDownIndices2, f_zero);
    } else {
        Bootstrap_RangeCheck_PatersonStockmeyer(range_check_res, result, rangeCheckIndices, p, ring_dim,
                                                relin_keys, seal_context, bfv_secret_key, f_zero, gateEval, skip_first_odd,
                                                firstDegree, secondDegree);
    }
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    cout << "TOTAL TIME for rangecheck: " << total_online << endl;
    cout << "Noise after rangecheck: " << decryptor.invariant_noise_budget(range_check_res) << " bits\n";

    // decryptor.decrypt(range_check_res, pl);
    // batch_encoder.decode(pl, v);
    // cout << "Decrypt after rangeCheck: \n" << v << endl;

    ////////////////////////////////////////// SLOT TO COEFFICIENT /////////////////////////////////////////////////////

    time_start = chrono::high_resolution_clock::now();
    evaluator.mod_switch_to_next_inplace(range_check_res);
    cout << "Noise after range check mod switch to very last???: " << decryptor.invariant_noise_budget(range_check_res) << " bits\n";

    Ciphertext range_check_res_copy(range_check_res);

    Evaluator eval_coeff(seal_context_last);
    eval_coeff.rotate_columns_inplace(range_check_res_copy, gal_keys_coeff);
    for (int i = 0; i < sq_ct; i++) {
        eval_coeff.rotate_rows(range_check_res, sq_rt * i, gal_keys_coeff, ct_sqrt_list[i]);
        eval_coeff.transform_to_ntt_inplace(ct_sqrt_list[i]);
        eval_coeff.rotate_rows(range_check_res_copy, sq_rt * i, gal_keys_coeff, ct_sqrt_list[i+sq_ct]);
        eval_coeff.transform_to_ntt_inplace(ct_sqrt_list[i+sq_ct]);
    }

    // vector<Plaintext> U_plain_list(ring_dim);
    // for (int iter = 0; iter < sq_rt; iter++) {
    //     for (int j = 0; j < (int) ct_sqrt_list.size(); j++) {
    //         vector<uint64_t> U_tmp = readUtemp(j*sq_rt + iter);
    //         batch_encoder.encode(U_tmp, U_plain_list[iter * ct_sqrt_list.size() + j]);
    //         evaluator.transform_to_ntt_inplace(U_plain_list[iter * ct_sqrt_list.size() + j], ct_sqrt_list[j].parms_id());
    //     }
    // }
    time_end = chrono::high_resolution_clock::now();
    total_preprocess += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();

    time_start = chrono::high_resolution_clock::now();
    // Ciphertext coeff = slotToCoeff(seal_context, seal_context_last, ct_sqrt_list, U_plain_list, gal_keys_coeff, sq_rt, ring_dim);
    Ciphertext coeff = slotToCoeff_WOPrepreocess(seal_context, seal_context_last, ct_sqrt_list, gal_keys_coeff, sq_rt, ring_dim);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    cout << "TOTAL TIME for slotToCoeff: " << total_online << endl;

    // cout << "PLAINTEX OF SLOTTOCOEFF\n";
    // decryptor.decrypt(coeff, pl);
    // for (int i = 0; i < ring_dim; i++) {
    //   cout << pl[i] << ",";
    // }
    // cout << endl;

    ////////////////////////////////////////////////// KEY SWITCHING ///////////////////////////////////////////////////

    time_start = chrono::high_resolution_clock::now();

    while(seal_context.last_parms_id() != coeff.parms_id()) {
        evaluator.mod_switch_to_next_inplace(coeff);
    }
    cout << "Noise before key switch: " << decryptor.invariant_noise_budget(coeff) << " bits\n";

    Ciphertext copy_coeff = coeff;
    auto ct_in_iter = util::iter(copy_coeff);
    ct_in_iter += coeff.size() - 1;
    seal::util::set_zero_poly(ring_dim, 1, coeff.data(1)); // notice that the coeff_mod.size() is hardcoded to 1, thus this needs to be performed on the last level

    evaluator.switch_key_inplace(coeff, *ct_in_iter, static_cast<const KSwitchKeys &>(ksk), 0, my_pool);

    cout << "Noise before extraction: " << decryptor.invariant_noise_budget(coeff) << " bits\n";
    
    vector<regevCiphertext> lwe_ct_results = extractRLWECiphertextToLWECiphertext(coeff);
    time_end = chrono::high_resolution_clock::now();
    total_online += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();
    cout << "TOTAL TIME for Extraction: " << total_online << endl;
    cout << "TOTAL PREPROCESS TIME: " << total_preprocess << endl;
    cout << "TOTAL ONLINE TIME: " << total_online << endl;

    return lwe_ct_results;
}


inline void multiply_power_of_X(EncryptionParameters& enc_param, const Ciphertext &encrypted, Ciphertext &destination, uint32_t index) {

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


// for a tree with m leaf node, m >> stepSize, we first expand it to a subtree with m / stepSize leaf node
// (i.e., this subtree is the top of the whole tree)
// and then for each leaf node in this subtree, expand it into a small subtree with stepSize leaf node
inline vector<Ciphertext> expand(const SEALContext& context, EncryptionParameters& enc_param,
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

    Plaintext pp;
    // for (int c = 0; c< (int) temp.size(); c++) {
    //   decryptor.decrypt(temp[c], pp);
    //   for (int i = 0; i < (int) m; i++) {
    //     cout << pp.data()[i] << " ";
    //   }
    //   cout << endl;
    // }
    // evaluator.add

    // Last step of the loop
    vector<Ciphertext> newtemp(temp.size() << 1);
    int index_raw = (m << 1) - (1 << (logm - 1));
    int index = (index_raw * galois_elts[logm - 1]) % (m << 1);

    for (uint32_t a = 0; a < temp.size(); a++) {
        if (a >= (m - (1 << (logm - 1)))) { // corner case.
            evaluator.multiply_plain(temp[a], two, newtemp[a]); // plain multiplication by 2.
        } else {
            evaluator.apply_galois(temp[a], galois_elts[logm - 1], galkey, tempctxt_rotated);

	    // cout << "Right before ??? " << a << endl;
	    // decryptor.decrypt(tempctxt_rotated, pp);
	    // for (int i = 0; i < (int) m; i++) {
	    //   cout << pp.data()[i] << " ";
	    // }
	    // cout << endl;
	    
            evaluator.add(temp[a], tempctxt_rotated, newtemp[a]);
            multiply_power_of_X(enc_param, temp[a], tempctxt_shifted, index_raw);
            multiply_power_of_X(enc_param, tempctxt_rotated, tempctxt_rotatedshifted, index);
            evaluator.add(tempctxt_shifted, tempctxt_rotatedshifted, newtemp[a + temp.size()]);
        }
    }

    vector<Ciphertext>::const_iterator first = newtemp.begin();
    vector<Ciphertext>::const_iterator last = newtemp.begin() + m;
    vector<Ciphertext> newVec(first, last);

    // cout << "********************************************\n";
    // for (int c = 0; c< (int) newVec.size(); c++) {
    //   decryptor.decrypt(newVec[c], pp);
    //   for (int i = 0; i < (int) m; i++) {
    // 	cout << pp.data()[i] << " ";
    //   }
    //   cout << endl;
    // }

    return newVec;
}

