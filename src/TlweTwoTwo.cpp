#include <iostream>
#include <math.h>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <tfhe/lwe-functions.h>
#include <tfhe/numeric_functions.h>
#include <tfhe/tlwe_functions.h>
#include <random>
#define MOD_MULT 64
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
bool TLweTwoPartyWithSmudge(int msg, TLweParams *params, TLweKey *key, double bound)
{
	std::cout << "N: 2^" << log2(params->N) << " Error bound: " << bound << "\n==================" << std::endl;
    int MSIZE = pow(2, ceil(log2(msg)));
    std::cout << "Plaintext: " << msg << "      MSIZE: " << MSIZE << std::endl;
    TLweSample *ciphertext_dir = new TLweSample(params);
    tLweSymEncryptT(ciphertext_dir, modSwitchToTorus32(msg,MSIZE), .001, key);
    TLweSample *ciphertext_thr = new TLweSample(params);
    tLweSymEncryptT(ciphertext_thr, modSwitchToTorus32(msg,MSIZE), .001, key);
    
    int dmsg_dir = modSwitchFromTorus32(tLweSymDecryptT(ciphertext_dir, key, MSIZE), MSIZE);
    std::cout << "After direct decryption: " << dmsg_dir << std::endl;

    TLweKey *secret1 = new TLweKey(params);
    tLweKeyGen(secret1);
    TorusPolynomial *intermediate_msg = new TorusPolynomial(params->N);
    tLwePhase(intermediate_msg, ciphertext_thr, secret1);
    for(int i = 0; i < params->N; i++){
    	ciphertext_thr->b->coefsT[i] = gaussian32(intermediate_msg->coefsT[i], bound);
    	// if(i<10) std::cout<<(ciphertext_thr->b->coefsT[i]-intermediate_msg->coefsT[i])<<std::endl;
    }
    for (int i = 0; i < params->k; i++){
        for (int j = 0; j < params->N; j++){
            key->key[i].coefs[j] -= secret1->key[i].coefs[j];
        }
    }
    int dmsg_thr = modSwitchFromTorus32(tLweSymDecryptT(ciphertext_thr, key, MSIZE), MSIZE);
    std::cout << "After 2-2 Threshold Decryption: " << dmsg_thr << "\n--------" << std::endl;
    return (msg == dmsg_thr);

    
}

int main()
{

    FILE *plaintext = fopen("test/plain22.txt", "r");
    int32_t msg;
    double b;
    fscanf(plaintext, "%d", &msg);
    fclose(plaintext);
    std::cout << msg << "\n";
    for (int i = 1024; i < 65536; i = i << 1){
    // for (int i = 1048576; i < 16777216; i = i << 1){
        TLweParams *params2 = new TLweParams(1024 * i, 2, 0.01, 0.2);
        TLweKey *key2 = new TLweKey(params2);
        tLweKeyGen(key2);

        std::cout << "\n\n\nTLWE\n\n\n\n";
        if (!TLweTwoPartyWithSmudge(msg, params2, key2, 0.00)){
            std::cout << "Threshold Decryption with zero smudging failed. Jumping to next power of 2.\n";
            delete_TLweKey(key2);
            delete_TLweParams(params2);
            continue;
        }
        b = 1.0;
        while (b > 1e-2){
            b /= 2;
            TLweTwoPartyWithSmudge(msg, params2, key2, b);
        }

        delete_TLweKey(key2);
        delete_TLweParams(params2);
    }

    return 0;
}