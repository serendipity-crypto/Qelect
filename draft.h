    // vector<uint64_t> test_v(ring_dim);
    // for (int i = 0; i < ring_dim; i++) {
    //     test_v[i] = i;
    // }
    // Plaintext test_p;
    // Ciphertext test_c;
    // batch_encoder.encode(test_v, test_p);
    // encryptor.encrypt(test_p, test_c);

    // for (int i = 0; i < ring_dim; i++) {
    //     cout << test_c.data(0)[i] << " ";
    // }
    // cout << endl;
    // for (int i = 0; i < ring_dim; i++) {
    //     cout << test_c.data(1)[i] << " ";
    // }
    // cout << endl;

    // cout << endl;
    // decryptor.decrypt(test_c, test_p);
    // for (int i = 0; i < ring_dim; i++) {
    //     cout << test_p.data()[i] << " ";
    // }
    // cout << endl;
    // // batch_encoder.decode(test_p, test_v);
    // // cout << "\n\nDecrypted: " << test_v << endl << endl;

    // inverse_ntt_negacyclic_harvey(bfv_secret_key.data().data(), seal_context.key_context_data()->small_ntt_tables()[0]);
	// for (int i = 0; i < ring_dim; i++) {
	// 	cout << bfv_secret_key.data()[i] << " ";
	// }
	// cout << endl;
	// seal::util::RNSIter new_key_rns(bfv_secret_key.data().data(), ring_dim);
	// ntt_negacyclic_harvey(new_key_rns, coeff_modulus.size(), seal_context.key_context_data()->small_ntt_tables());











    // Plaintext blind_rot_base_p;
    // blind_rot_base_p.resize(ring_dim);
    // blind_rot_base_p.parms_id() = parms_id_zero;
    // for (int j = 0; j < ring_dim; j++) {
    //     blind_rot_base_p.data()[j] = 0;
    // }
    // blind_rot_base_p.data()[ring_dim-1] = p-1; // encode -X^{-1}, notice that X^{2N} = -1
    // evaluator.multiply_plain_inplace(final_perm_vec, blind_rot_base_p); // place the second element to the constant term












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
    // cout << "** [TIME] ** Expansion time: " << 
    //       chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() << endl;
    // total_time += chrono::duration_cast<chrono::microseconds>(time_end - time_start).count();


    // // multiply the selection binary ciphertexts with the tokens, and sum them up
    // time_start = chrono::high_resolution_clock::now();

    // vector<Ciphertext> result_tmp(numcores);

    // for (int i = 0; i < (int) expanded_leaf.size(); i++) {
    //     for (int j = 0; j < (int) expanded_leaf[0].size(); j++) {
    //         evaluator.mod_switch_to_next_inplace(expanded_leaf[i][j]);
    //         evaluator.transform_to_ntt_inplace(expanded_leaf[i][j]);
    //     }
    // }
    // evaluator.transform_to_ntt_inplace(tokens, expanded_leaf[0][0].parms_id());

    // int batch_share = group_size/numcores;
    // NTL_EXEC_RANGE(numcores, first, last);
    // for (int t = first; t < last; t++) {
    //     for (int i = t * batch_share; i < (t+1) * batch_share; i++) {
    //         if (i % batch_share == 0) {
    //             evaluator.multiply_plain(expanded_leaf[t][0], tokens, result_tmp[t]);
    //         } else {
    //             Ciphertext tmp;
    //             evaluator.multiply_plain(expanded_leaf[t][i % batch_share], tokens, tmp);
    //             evaluator.add_inplace(result_tmp[t], tmp);
    //         }
    //     }
    // }
    // NTL_EXEC_RANGE_END;
