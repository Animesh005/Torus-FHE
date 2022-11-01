#include <iostream>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <tfhe/lwe-functions.h>
#include <tfhe/numeric_functions.h>
#include <tfhe/tlwe_functions.h>
#include <random>
#define MSIZE 2

int n_smudge = 16;

int NSmudges(int n)
{
    return n_smudge;
}

void RandomSmudge(TorusPolynomial *b, double bound, int N, int r)
{
    int cnt = 0;
    int idx = 0;
    
    while (cnt < r){
        int off = rand() % N;
        if (off == 0)
            continue;

        idx = (idx + off) % N;
        b->coefsT[idx] += gaussian32(0, bound);
        cnt++;
    }
}


bool TLweNPartyWithSmudge(int msg, TLweParams *params, TLweKey *key, double bound, int n)
{
    std::cerr << "Error bound: " << bound << std::endl;
    std::cerr << "Plaintext: " << msg << std::endl;

    TLweSample **ciphertext = new TLweSample*[32];
    for (int i = 0; i < 32; i++){
        ciphertext[i] = new_TLweSample(params);
        tLweSymEncryptT(ciphertext[i], modSwitchToTorus32((msg >> i) & 1, MSIZE), 0.001, key);
    }

    int dmsg = 0;
    for (int i = 0; i < 32; i++){
        int bit = (tLweSymDecryptT(ciphertext[i], key, 2) == 0) ? 0 : 1;
        dmsg += (bit << i);
    }

    std::cerr << "Sanity check " << msg << " " << dmsg << std::endl;

    TLweKey *key2;
    TorusPolynomial **result;

    for (int k = 0; k < n - 1; k++){
        key2 = new_TLweKey(params);
        tLweKeyGen(key2);

        result = new TorusPolynomial*[32];

        for (int i = 0; i < 32; i++){
            result[i] = new_TorusPolynomial(params->N);
            tLwePhase(result[i], ciphertext[i], key2);
        }

        for (int i = 0; i < 32; i++){
            for (int j = 0; j < params->N; j++)
                ciphertext[i]->b->coefsT[j] = result[i]->coefsT[j];

            RandomSmudge(ciphertext[i]->b, bound, params->N, NSmudges(params->N));
        }

        std::cerr << "Decryption number: " << k << std::endl;

        for (int i = 0; i < params->k; i++){
            for (int j = 0; j < params->N; j++){
                key->key[i].coefs[j] -= key2->key[i].coefs[j];
            }
        }

        delete_TLweKey(key2);
        for (int i = 0; i < 32; i++)
            delete_TorusPolynomial(result[i]);
    }

    std::cerr << "Decryption number: " << n - 1 << std::endl;
    dmsg = 0;
    for (int i = 0; i < 32; i++){
        int bit = (tLweSymDecryptT(ciphertext[i], key, 2) == 0) ? 0 : 1;
        dmsg += (bit << i);
    }

    std::cerr << "Decrypted: " << dmsg << "\n\n";

    for (int i = 0; i < 32; i++){
        delete_TLweSample(ciphertext[i]);
    }

    return (msg == dmsg);
}


int main(int argc, char *argv[])
{
    n_smudge = atoi(argv[1]);

    FILE *plaintext = fopen("test/plain22.txt", "r");
    int32_t msg;
    fscanf(plaintext, "%d", &msg);
    fclose(plaintext);
    

    TLweParams *params = new_TLweParams(1024, 2, 0.01, 0.2);
    TLweKey *key;

    for (int p = 2; p <= 20; p++){
        for (double bound = 0.75; bound > 1e-6; bound /= 2){
            key = new_TLweKey(params);
            tLweKeyGen(key);
            srand(42);
            if (TLweNPartyWithSmudge(msg, params, key, bound, p)){
                std::cout << p << ", "<< bound << std::endl;
                break;
            }
        }
    }

    delete_TLweParams(params);

    return 0;
}