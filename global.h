#pragma once
#include "seal/seal.h"
using namespace seal;

size_t poly_modulus_degree_glb = 32768;

// (root, ring_dim) --> root^(2*ring_dim) % 65527 = 1
// (4, 8), (2, 16), (255, 32), (141, 128), (431, 512), (21, 2048), (15, 8192), (3, 32768)
int primitive_root = 3;

// plaintext_to_ciphertext_prime_map:
// equal and below 2^9 -> 65537, 2^12 -> (20bit prime) 786433 = 512*512*3+1, 2^15 -> (23bit prime) 5308417 = (2^8*9)^2+1
int prime_p = 786433;

int expand_core_share = 8192; // 32GB mem can only tolerate 8192 ciphertexts being generated at the same time

uint64_t loading_time;
uint64_t U_time;
uint64_t U_time_multi_core;
