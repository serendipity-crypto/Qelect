#include <NTL/BasicThreadPool.h>
#include <NTL/ZZ.h>
#include <thread>
#include "seal/util/polyarithsmallmod.h"
#include "seal/seal.h"

using namespace seal::util;
using namespace std;
using namespace seal;


inline
long power_seal(long x, long y, long m)
{
    if (y == 0)
        return 1;
    long p = power_seal(x, y / 2, m) % m;
    p = (p * p) % m;
 
    return (y % 2 == 0) ? p : (x * p) % m;
}

inline
long modInverse_seal(long a, long m)
{
    return power_seal(a, m - 2, m);
}


vector<vector<uint64_t>> generateLowerTri_vandermonde(vector<vector<uint64_t>> vander, const int q = 65537) {
    int size = (int) vander.size(); 
    vector<vector<uint64_t>> lower(size, vector<uint64_t>(size));

    // x_j = vander[j][1];
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i < j) {
                lower[i][j] = 0;
            } else if (i == 0 && j == 0) {
                lower[i][j] = 1;
            } else if (i == j){
                long tmp = 1;
                for (int k = 0; k <= i; k++) {
                    if (k != j) {
                        long ttmp = vander[j][1] - vander[k][1];
                        while (ttmp < 0) ttmp = q + ttmp;
                        ttmp = modInverse_seal(ttmp, q); // 1 / (x_j-x_k)
                        tmp = (tmp * ttmp) % q;
                    }
                }
                lower[i][j] = (uint64_t) (tmp % q);
            } else {
                long tmp = (long) lower[i-1][j];
                long ttmp = vander[j][1] - vander[i][1];
                while (ttmp < 0) ttmp = q + ttmp;
                ttmp = modInverse_seal(ttmp, q);
                tmp = (tmp * ttmp) % q;
                lower[i][j] = (uint64_t) (tmp % q);
            }
        }
    }

    return lower;
}


vector<vector<uint64_t>> generateUpperTri_vandermonde(vector<vector<uint64_t>> vander, const int q = 65537) {
    int size = (int) vander.size(); 
    vector<vector<uint64_t>> upper(size, vector<uint64_t>(size));

    // x_j = vander[j][1];
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i == j) {
                upper[i][j] = 1;
            } else if (j == 0) {
                upper[i][j] = 0;
            } else {
                long tmp = (i-1 < 0) ? - upper[i][j-1] * vander[j-1][1] : upper[i-1][j-1] - upper[i][j-1] * vander[j-1][1];
                while (tmp < 0) tmp = q + tmp;
                upper[i][j] = (uint64_t) (tmp % q);
            }
        }
    }

    return upper;
}

// compute C:= A x B, where A, B, C are matrices
vector<vector<uint64_t>> matrixMultiplication(vector<vector<uint64_t>> a, vector<vector<uint64_t>> b, const int q = 65537) {
    int size = (int) a.size(); 
    vector<vector<uint64_t>> result(size, vector<uint64_t>(size));

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            long temp = 0;
            for (int k = 0; k < (int) result.size(); k++) {
                temp = (temp + a[i][k] * b[k][j]) % q;
                temp = temp < 0 ? temp + q : temp;
            }
            result[i][j] = (uint64_t) temp;
        }
    }

    return result;
}

vector<vector<uint64_t>> generateInverse_vander(vector<vector<uint64_t>> vander, const int q = 65537) {
    vector<vector<uint64_t>> upper = generateUpperTri_vandermonde(vander);
    cout << "       upper done...\n";
    vector<vector<uint64_t>> lower = generateLowerTri_vandermonde(vander);
    cout << "       lower done...\n";

    vector<vector<uint64_t>> result = matrixMultiplication(upper, lower);

    return result;
}