#include <iostream>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <tfhe/lwe-functions.h>
#include <tfhe/numeric_functions.h>
#include <tfhe/tlwe_functions.h>
#include <random>

void printTorusPolynomial(TorusPolynomial *p){
    std::cout << "[ ";

    for (int i = 0; i < p->N; i++){
        std::cout << p->coefsT[i] << " ";
    }

    std::cout << "]\n";
}

void printIntPolynomial(IntPolynomial *p){
    std::cout << "[ ";

    for (int i = 0; i < p->N; i++){
        std::cout << p->coefs[i] << " ";
    }

    std::cout << "]\n";
}

void nonFFTmul(TorusPolynomial *ans, IntPolynomial *S, TorusPolynomial *A, int mod)
{
    long long *temp = new long long[S->N + A->N];
    for (int i = 0; i < S->N + A->N; i++)
        temp[i] = 0;

    for (int i = 0; i < S->N; i++){
        for (int j = 0; j < A->N; j++){
            temp[i + j] += ((int64_t)S->coefs[i] * (int64_t)A->coefsT[j]); 
        }
    }

    for (int i = mod; i < S->N + A->N; i++){
        temp[i % mod] -= temp[i];
    }

    for (int i = 0; i < mod; i++){
        ans->coefsT[i] = (temp[i] % ((long long)549755809793));
    }
}

int main()
{
    TorusPolynomial *A = new_TorusPolynomial(1024);
    IntPolynomial *S = new_IntPolynomial(1024);
    TorusPolynomial *ans = new_TorusPolynomial(1024);
    
    std::default_random_engine generator;
    // generator.seed(142);
    std::uniform_int_distribution<int> distribution(0,9);
    // distribution(generator);
    
    for (int i = 0; i < 1024; i++){
        A->coefsT[i] = gaussian32(0, 3.0);
        S->coefs[i] = distribution(generator);
        ans->coefsT[i] = 0;
    }


    torusPolynomialAddMulR(ans, S, A);

    std::cout << "A: ";
    printTorusPolynomial(A);

    std::cout << "S: ";
    printIntPolynomial(S);

    std::cout << "ans: ";
    printTorusPolynomial(ans);

    TorusPolynomial *ans2 = new_TorusPolynomial(1024);
    nonFFTmul(ans2, S, A, 1024);

    std::cout << "ans2: ";
    printTorusPolynomial(ans2);

    for (int i = 0; i < 1024; i++){
        ans->coefsT[i] = (ans->coefsT[i] & 2147483647);
        ans2->coefsT[i] = (ans2->coefsT[i] & 2147483647);
        
        if (ans->coefsT[i] != ans2->coefsT[i])
            std::cout << "Mismatch: " << i << " " << ans->coefsT[i] - ans2->coefsT[i] << std::endl;
    }

    return 0;
}