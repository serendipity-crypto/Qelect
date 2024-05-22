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











