#include <iostream>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <tfhe/lwe-functions.h>
#include <tfhe/numeric_functions.h>
#include <tfhe/tlwe_functions.h>
#include <random>
#include "share.hpp"

void onesComplemet(LweSample *result, const LweSample *input1, const int nbits, const TFheGateBootstrappingCloudKeySet *bk)
{
    for (int i = 0; i < nbits; i++){
        bootsNOT(&result[i], &input1[i], bk);
        //bootsXOR(&result[i], &input1[i], &input2[i], bk);
    }
}


void Compute(LweSample *result, const LweSample *input1, const LweSample *input2, const int nbits, const TFheGateBootstrappingCloudKeySet *bk)
{
    for (int i = 0; i < nbits; i++){
        bootsXOR(&result[i], &input1[i], &input2[i], bk);
    }
}

bool TOutOfN(int msg, int msg2, TFheGateBootstrappingParameterSet *params, TFheGateBootstrappingSecretKeySet *key, const TFheGateBootstrappingCloudKeySet *cloud_key, int t, int n)
{
    std::cout << "Plaintext: " << msg << std::endl;
    std::cout << "Plaintext2: " << msg2 << std::endl;
    int resultant_msg = msg^msg2;
    std::cout << "Resultant_Plaintext: " << resultant_msg<< std::endl;


    LweSample *ciphertext = new_gate_bootstrapping_ciphertext_array(32, params);
    for (int i = 0; i < 32; i++)
        bootsSymEncrypt(&ciphertext[i], (msg >> i) & 1, key);
    int sz = key->lwe_key->params->n;

    LweSample *ciphertext2 = new_gate_bootstrapping_ciphertext_array(32, params);
    for (int i = 0; i < 32; i++)
        bootsSymEncrypt(&ciphertext2[i], (msg2 >> i) & 1, key);

    LweSample *result = new_gate_bootstrapping_ciphertext_array(32, params);
    Compute(result, ciphertext, ciphertext2, 32, cloud_key);
    


    ublas::vector<INT> sk(sz);
    for (int i = 0; i < sz; i++){
        sk(i) = key->lwe_key->key[i];
    }

    auto shares = gen_shares(sk, t, n, 3);
    int dmsg = 0;
    for (int i = 0; i < 32; i++){
        ublas::vector<INT> a(sz);
        for (int j = 0; j < sz; j++){
            a(j) = result[i].a[j];
        }
        auto pdts = apply_product(shares, a);
        Torus32 axs = recostruct_combination(pdts, t);
        Torus32 phase = result[i].b - axs;
        int bit = (phase > 0) ? 1 : 0;
        dmsg += (bit << i);
    }
    delete_gate_bootstrapping_ciphertext_array(32, ciphertext);
    delete_gate_bootstrapping_ciphertext_array(32, ciphertext2);
    delete_gate_bootstrapping_ciphertext_array(32, result);

    std::cout << "Direct decryption: " << resultant_msg << "\nThreshold Decryption: " << dmsg << std::endl;
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

    FILE *plaintext2 = fopen("test/plain23.txt", "r");
    int32_t msg2;
    fscanf(plaintext2, "%d", &msg2);
    fclose(plaintext2);
    
    
    FILE *paramFile = fopen("test/secret.params", "rb");
    auto params = new_tfheGateBootstrappingParameterSet_fromFile(paramFile);
    fclose(paramFile);

    FILE *cloudKeyFile = fopen("test/cloud.key", "rb");
    auto cloud_key_ = new_tfheGateBootstrappingCloudKeySet_fromFile(cloudKeyFile);
    fclose(cloudKeyFile);    

    TOutOfN(msg, msg2, params, key, cloud_key_, 2, 3);
    return 0;
}
