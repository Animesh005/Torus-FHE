#include <iostream>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <tfhe/lwe-functions.h>
#include <tfhe/numeric_functions.h>
#include <tfhe/tlwe_functions.h>
#include <random>
#include "share.hpp"

bool TOutOfN(int msg, TFheGateBootstrappingParameterSet *params, TFheGateBootstrappingSecretKeySet *key, int t, int n)
{
    std::cout << "Plaintext: " << msg << std::endl;

    LweSample *ciphertext = new_gate_bootstrapping_ciphertext_array(32, params);
    for (int i = 0; i < 32; i++)
        bootsSymEncrypt(&ciphertext[i], (msg >> i) & 1, key);
    int sz = key->lwe_key->params->n;

    ublas::vector<INT> sk(sz);
    for (int i = 0; i < sz; i++){
        sk(i) = key->lwe_key->key[i];
    }

    auto shares = gen_shares(sk, t, n, 3);
    int dmsg = 0;
    for (int i = 0; i < 32; i++){
        ublas::vector<INT> a(sz);
        for (int j = 0; j < sz; j++){
            a(j) = ciphertext[i].a[j];
        }
        auto pdts = apply_product(shares, a);
        Torus32 axs = recostruct_combination(pdts, t);
        Torus32 phase = ciphertext[i].b - axs;
        int bit = (phase > 0) ? 1 : 0;
        dmsg += (bit << i);
    }
    delete_gate_bootstrapping_ciphertext_array(32, ciphertext);

    std::cout << "Direct decryption: " << msg << "\nThreshold Decryption: " << dmsg << std::endl;
    return (msg == dmsg);
}

int main()
{
    // Read from file as binary and encrypt
    FILE *skFile = fopen("test/secret.key", "rb");
    auto key = new_tfheGateBootstrappingSecretKeySet_fromFile(skFile);
    fclose(skFile);

    FILE *plaintext = fopen("test/plain22.txt", "r");
    int32_t msg;
    fscanf(plaintext, "%d", &msg);
    fclose(plaintext);
    
    
    FILE *paramFile = fopen("test/secret.params", "rb");
    auto params = new_tfheGateBootstrappingParameterSet_fromFile(paramFile);
    fclose(paramFile);    

    TOutOfN(msg, params, key, 2, 3);
    return 0;
}
