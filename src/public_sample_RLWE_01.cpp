#include <iostream>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <tfhe/lwe-functions.h>
#include <tfhe/numeric_functions.h>
#include <tfhe/tfhe_garbage_collector.h>
#include <random>
#include <time.h>

typedef struct myParameterSet{
	TLweParams *tlweparams;
	TFheGateBootstrappingParameterSet *bootparams;
}myParameterSet;
myParameterSet *initialize_gate_bootstrapping_params() {
    static const int32_t N = 1024;
    static const int32_t k = 1;
    static const int32_t n = 1024;
    static const int32_t bk_l = 3;
    static const int32_t bk_Bgbit = 7;
    static const int32_t ks_basebit = 2;
    static const int32_t ks_length = 8;
    static const double ks_stdev = pow(2.,-15); //standard deviation
    static const double bk_stdev = pow(2.,-25);; //standard deviation
    static const double max_stdev = 0.012467; //max standard deviation for a 1/4 msg space

    LweParams *params_in = new_LweParams(n, ks_stdev, max_stdev);
    TLweParams *params_accum = new_TLweParams(N, k, bk_stdev, max_stdev);
    TGswParams *params_bk = new_TGswParams(bk_l, bk_Bgbit, params_accum);

    TfheGarbageCollector::register_param(params_in);
    TfheGarbageCollector::register_param(params_accum);
    TfheGarbageCollector::register_param(params_bk);

    auto bootparams = new TFheGateBootstrappingParameterSet(ks_length, ks_basebit, params_in, params_bk);
    myParameterSet *params = new myParameterSet();
    params->tlweparams = params_accum;
    params->bootparams = bootparams;
    return params;
}

void convert_rlwe_to_lwe(LweSample *c, TLweSample *rc, const LweParams *lweparams, int pos){
	int n = lweparams->n;
	c->b = rc->b->coefsT[pos];
	for(int i = 0, j = pos; j >= 0; j--, i++){
		c->a[i] = rc->a[0].coefsT[j];
	}
	for(int i = pos+1, j = n-1; j >= pos+1; j--, i++){
		c->a[i] = -(rc->a[0].coefsT[j]);
	}
}

void convert_lwe_to_rlwe(TLweSample *rc, LweSample *c, const TLweParams *tlweparams){
	int N = tlweparams->N;
	rc->b->coefsT[0] = c->b;
	rc->a[0].coefsT[0] = c->a[0];
	for(int i = 1; i < N; i++){
		rc->a[0].coefsT[i] = -(c->a[N-i]);
	}
}

void add_msg(TLweSample *result, int msg){
	if(msg == 0){
		result->b->coefsT[0] = result->b->coefsT[0] + modSwitchToTorus32(0,2);
	}
	else{
		result->b->coefsT[0] = result->b->coefsT[0] + modSwitchToTorus32(1,2);
	}
}

void generate_sample(TLweSample *sample, TLweSample *input, LweSample *intermediate, const LweParams *lweparams, const TLweParams *tlweparams, const TFheGateBootstrappingCloudKeySet *cloud_key, int bit){
	LweSample *input_lwe = new_LweSample(lweparams);
	convert_rlwe_to_lwe(input_lwe, input, lweparams, 0);
	bootsXOR(intermediate, input_lwe, input_lwe, cloud_key);
	convert_lwe_to_rlwe(sample, intermediate, tlweparams);
	add_msg(sample, bit);
}

int main(){
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dist(0, 1);

	auto myparams = initialize_gate_bootstrapping_params();
	const TLweParams *tlweparams = myparams->tlweparams;
	auto params = myparams->bootparams;
	const LweParams *lweparams = params->in_out_params;

	TFheGateBootstrappingSecretKeySet *key = new_random_gate_bootstrapping_secret_keyset(params);
	auto cloud_key = &key->cloud;

	/* choose a random bit */
	int sel = dist(gen);

	/* TLWE key from LWE key */
	TLweKey *tlwekey = new_TLweKey(tlweparams);
	for(int i = 0; i < tlweparams->N; i++){
		tlwekey->key[0].coefs[i] = key->lwe_key->key[i];
	}

	/* Encryption : Encrypt message 0 if sel = 0, else encrypt message 1 */
	TLweSample *ciphertext = new_TLweSample(tlweparams);
	if(sel == 0){
		tLweSymEncryptT(ciphertext, -modSwitchToTorus32(1, 8), 0.001, tlwekey);
	}
	else{
		tLweSymEncryptT(ciphertext, modSwitchToTorus32(1, 8), 0.001, tlwekey);
	}

	std::cout << "pass encryption of: " << sel << "\n\n"; 
	
	/* Generate fresh samples from two distributions, first from the encryptions of 0 and second from the encryptions of 1 */
	TLweSample *fresh_sample0 = new_TLweSample(tlweparams);
	TLweSample *fresh_sample1 = new_TLweSample(tlweparams);

	LweSample *intermediate = new_LweSample(lweparams);

	/* verify correctness of sampling algorithm */
	std::cout << "Generating sample of encryption of 0\n";
	generate_sample(fresh_sample0, ciphertext, intermediate, lweparams, tlweparams, cloud_key, 0);
	int dmsg0 = tLweSymDecryptT(fresh_sample0, tlwekey, 8) > 0 ? 1 : 0;
	std::cout << "expected: 0, generated: " << dmsg0 << "\n\n";

	std::cout << "Generating sample of encryption of 1\n";
	generate_sample(fresh_sample1, ciphertext, intermediate, lweparams, tlweparams, cloud_key, 1);
	int dmsg1 = tLweSymDecryptT(fresh_sample1, tlwekey, 8) > 0 ? 1 : 0;
	std::cout << "expected: 1, generated: " << dmsg1 << "\n\n";
	
	return 0;
}