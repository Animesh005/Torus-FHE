#include <iostream>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <tfhe/lwe-functions.h>
#include <tfhe/numeric_functions.h>
#include <tfhe/tlwe_functions.h>
#include <random>
#include <time.h>
#include <cstdint>
#include <x86intrin.h>
#include <bits/stdc++.h>
#include "threshold_decryption_functions.hpp"


extern std::map<std::pair<int, int>, int> ncr_cacheT;	/* Stores <<n, r>: C(n, r)> */
extern std::map<std::pair<int, int>, TLweKey*> shared_key_repo;	/* Stores <<party_id, group_id>: key_share> */



#define MSIZE 2
// extern std::map<std::pair<int, int>, TLweKey*> shared_key_repo;

int main(int argc, char *argv[])
{
    if(argc < 3){
		std::cerr << "Please provide values of t, p in the command line for t-out-of-p threshold decryption.\n";
		return 0;
	}

	int t = atoi(argv[1]);
	int p = atoi(argv[2]);

    FILE *plaintext = fopen("test/plain22.txt", "r");
    int32_t msg;
    fscanf(plaintext, "%d", &msg);
    fclose(plaintext);

    /* Set Up */
	TLweParams *params = new_TLweParams(1024, 1, 0.01, 0.2);
	TLweKey *key = new_TLweKey(params);
	tLweKeyGen(key);

	/* Encryption */
	TLweSample *ciphertext = new_TLweSample(params);
	TorusPolynomial* mu = new_TorusPolynomial(params->N);
	for(int i = 0; i < params->N; i++){
		mu->coefsT[i] = 0;
	}
	for(int i = 0; i < 32; i++){
		mu->coefsT[i] += modSwitchToTorus32((msg >> i) & 1, MSIZE);
	}
	tLweSymEncrypt(ciphertext, mu, 0.001, key);

    std::cout << "#pragma once\n#include \"tfhe/common.h\"\n#include \"tfhe/rlwe.h\"\n";
    std::cout << "TLweParams params = {\n";
    std::cout << "\t.N = " << params->N << ",\n";
    std::cout << "\t.k = " << params->k << ",\n";
    std::cout << "\t.alpha_max = " << params->alpha_max << ",\n";
    std::cout << "\t.alpha_min = " << params->alpha_min << ",\n";
    std::cout << "\t.extracted_lweparams = {\n"; 
    std::cout << "\t\t.alpha_max = " << params->extracted_lweparams.alpha_max << ",\n";
    std::cout << "\t\t.alpha_min = " << params->extracted_lweparams.alpha_min << ",\n";
    std::cout << "\t\t.n = " << params->extracted_lweparams.n << "\n";
    std::cout << "\t}" << std::endl;
    std::cout << "};" << std::endl;

    int n = params->N;
    
    std::cout << "const Torus32 __adata[] = {\n";
    for (int i = 0; i < n; i++){
        if (i % 10 == 0){
            std::cout << "\t";
        }
        std::cout << ciphertext->a->coefsT[i] << ",";
        if (i % 10 == 9){
            std::cout << "\n";
        }else{
            std::cout << " ";
        }
    }
    std::cout << "\n};" << std::endl;

    std::cout << "const Torus32 __bdata[] = {\n";
    for (int i = 0; i < n; i++){
        if (i % 10 == 0){
            std::cout << "\t";
        }
        std::cout << ciphertext->b->coefsT[i] << ",";
        if (i % 10 == 9){
            std::cout << "\n";
        }else{
            std::cout << " ";
        }
    }
    std::cout << "\n};" << std::endl;

    shareSecret(t, p, key, params);

    int ncr = ncrT(p, t);

    // int iii = 0;
    // for (auto it: shared_key_repo){
    //     std::cout << it.first.first << ", " << it.first.second << std::endl;
    //     iii++;
    // }
    // std::cout << ncr << " " << iii << std::endl;
    
    // Switching from 1-based to 0-based index
    std::cout << "const int32_t __keydata_coef[" << p << "][" << ncr << "][" << n << "] = {\n";
    for (int i = 0; i < p; i++){
        std::cout << "\t{\n";
        for (int j = 0; j < ncr; j++){
            std::cout << "\t\t{\n";
            if (shared_key_repo.find({i + 1, j + 1}) != shared_key_repo.end()){
                TLweKey *_key = shared_key_repo[{i + 1, j + 1}];
                for (int k = 0; k < n; k++){
                    if (k % 10 == 0){
                        std::cout << "\t\t\t";
                    }
                    std::cout << _key->key->coefs[k] << ",";
                    if (k % 10 == 9){
                        std::cout << "\n";
                    }else{
                        std::cout << " ";
                    }
                }
            }
            std::cout << "\t\t},\n";
        }
        std::cout << "\t},\n";

    }
    std::cout << "};" << std::endl;

    std::cout << "const Torus32 zero[] = {\n";
    for (int i = 0; i < n; i++){
        if (i % 10 == 0){
            std::cout << "\t";
        }
        std::cout << 0 << ",";
        if (i % 10 == 9){
            std::cout << "\n";
        }else{
            std::cout << " ";
        }
    }
    std::cout << "};" << std::endl;

    std::cout << "TLweSample result = {\n";
    std::cout << "\t.current_variance = " << ciphertext->current_variance << "\n";
    std::cout << "};" << std::endl;

    return 0;
}