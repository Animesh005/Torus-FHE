#include <iostream>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <tfhe/lwe-functions.h>
#include <tfhe/numeric_functions.h>
#include <tfhe/tlwe_functions.h>
#include <random>
#include <time.h>
#define MSIZE 2
#define MOD_MULT 64

/**
 * b - a.s a Z_q s 0, 1
 * b = a.s + e
 * s = s1 + s2
 * x = b - a.s1 + e1
 * y = x - a.s2
 */

int  count =0;
clock_t times =0;


bool TwoPartyWithSmudge(int msg, TFheGateBootstrappingParameterSet *params, TFheGateBootstrappingSecretKeySet *key, double bound)
{
    std::cout << "Error bound: " << bound << "\n==================" << std::endl;
    std::cout << "Plaintext: " << msg << std::endl;

    LweSample *ciphertext = new_gate_bootstrapping_ciphertext_array(32, params);
    for (int i = 0; i < 32; i++)
        bootsSymEncrypt(&ciphertext[i], (msg >> i) & 1, key);
    int sz = key->lwe_key->params->n;
    int32_t *old_key = new int32_t[sz];
    int32_t *new_key1 = new int32_t[sz];

    std:: cout<<"\n n ="<< key->lwe_key->params->n;

    clock_t start2,end2;
    start2= clock();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, 1);

    for (int i = 0; i < sz; i++){
        old_key[i] = key->lwe_key->key[i];
        new_key1[i] = dist(gen);
        // std:: cout<<"\n old key[]="<<old_key[i];
        // std:: cout<<"\n new key[]="<<new_key1[i];
        key->lwe_key->key[i] = new_key1[i];
    }

    // std:: cout << old_key[0]<<"\n";
    // std:: cout << old_key[1]<<"\n";
    // std:: cout << new_key1[0]<<"\n";
    // std:: cout << new_key1[1]<<"\n";

    // First partial decryption
    for (int i = 0; i < 32; i++){
        if (bound > 0)
            ciphertext[i].b = gaussian32(lwePhase(&ciphertext[i], key->lwe_key), bound);
        else
            ciphertext[i].b = lwePhase(&ciphertext[i], key->lwe_key);
        
    }
    std::cout << "\nFirst Partial Decryption done" << std::endl;
    // Second partial decryption
    for (int i = 0; i < sz; i++){
        key->lwe_key->key[i] = old_key[i] - new_key1[i];
    }

    int dmsg = 0;
    for (int i = 0; i < 32; i++){
        int bit = bootsSymDecrypt(&ciphertext[i], key);
        dmsg += (bit << i);
    }

    end2 = clock();
    times = times + (end2-start2);
    count = count + 1;
    printf("Two decrypt time : %.3f ms\n", (double(end2-start2)/CLOCKS_PER_SEC)*1000);
    std::cout << "Decrypted: " << dmsg << "\n\n";

    delete_gate_bootstrapping_ciphertext_array(32, ciphertext);

    return (msg == dmsg);
}

bool TLweTwoPartyWithSmudge(int msg, TLweParams *params, TLweKey *key, double bound)
{
    std::cout << "Error bound: " << bound << "\n==================" << std::endl;
    std::cout << "Plaintext: " << msg << std::endl;

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

    std::cout << "Sanity check " << msg << " " << dmsg << std::endl;

    auto key2 = new_TLweKey(params);
    tLweKeyGen(key2);
    std::cout<<"N = "<<key2->params->N<<"\n";
    std::cout<<"K = "<<key2->params->k<<"\n";
    std::cout <<"key2 = "<<key2;

    for (int i = 0; i < params->k; i++){
        for (int j = 0; j < 10; j++){
            std::cout << key->key[i].coefs[j] << " ";
        }
        std::cout << "\n";
    }

    TorusPolynomial **result = new TorusPolynomial*[32];

    // 1st decryption
    for (int i = 0; i < 32; i++){
        result[i] = new_TorusPolynomial(params->N);
        tLwePhase(result[i], ciphertext[i], key2);
        // tLweApproxPhase(result[i], result[i], MSIZE, params->N);
    }

    for (int i = 0; i < 32; i++){
        for (int j = 0; j < params->N; j++)
            ciphertext[i]->b->coefsT[j] = result[i]->coefsT[j] + dtot32(bound);
    }

    // 10 10 10 10 10 10 .... N
    // m

    std::cout << "First Decryption Done" << std::endl;

    for (int i = 0; i < params->k; i++){
        for (int j = 0; j < params->N; j++){
            //std::cout<<key->key[i].coefs[j]<<" "<<key2->key[i].coefs[j]<<"\n";
            key->key[i].coefs[j] -= key2->key[i].coefs[j];
           // std::cout<<key->key[i].coefs[j]<<" "<<key2->key[i].coefs[j]<<"\n\n";
        }
    }

    dmsg = 0;
    for (int i = 0; i < 32; i++){
        int bit = (tLweSymDecryptT(ciphertext[i], key, 2) == 0) ? 0 : 1;
        dmsg += (bit << i);
    }

    // Restore the key
    for (int i = 0; i < params->k; i++){
        for (int j = 0; j < params->N; j++){
            key->key[i].coefs[j] += key2->key[i].coefs[j];
        }
    }

    std::cout << "Decrypted: " << dmsg << "\n\n";

    for (int i = 0; i < 32; i++){
        delete_TorusPolynomial(result[i]);
        delete_TLweSample(ciphertext[i]);
    }

    return (msg == dmsg);
}


// bool TGswTwoPartyWithSmudge(int msg, TGswParams *params, TGswKey *key, double bound)
// {
//     TGswSample *ciphertext = new_TGswSample(params);
//     tGswSymEncryptInt(ciphertext, msg, 0.42, key);

//     auto key1 = new_TGswKey(params);

//     tGswSymDecryptInt()

// }

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

    std::cout << "LWE\n\n\n\n";
    TwoPartyWithSmudge(msg, params, key, 0);
    double b = 1.0;
    while (b > 1e-2){
        b /= 2;
        TwoPartyWithSmudge(msg, params, key, b);
    }

    printf("Avg TwoTwo decrypt time : %.3f ms\n", ((double(times)/count)/CLOCKS_PER_SEC)*1000);


    // std::cout << "TLWE\n\n\n\n";

    // TLweParams *params2;
    // TLweKey *key2;

    // for (int i = 128; i <= 1024; i = i << 1){
    //     params2 = new_TLweParams(1024 * i, 2, 0.01, 0.2);
    //     key2 = new_TLweKey(params2);
    //     tLweKeyGen(key2);

    //     std::cout << "\n\nN = " << 1024 * i << std::endl; 
    //     if (!TLweTwoPartyWithSmudge(msg, params2, key2, 0)){
    //         std::cout << 1024 * i << " " << 0 << std::endl;
    //         continue;
    //     }
    //     b = 1.0;
    //     while (b > 1e-2){
    //         b /= 2;
    //         if (TLweTwoPartyWithSmudge(msg, params2, key2, b)){
    //             std::cout << 1024 * i << " " << b << std::endl;
    //             break;
    //         }
    //     }

    //     delete_TLweKey(key2);
    //     delete_TLweParams(params2);
    // }

    // delete_gate_bootstrapping_secret_keyset(key);
    // delete_gate_bootstrapping_parameters(params);

    return 0;
}