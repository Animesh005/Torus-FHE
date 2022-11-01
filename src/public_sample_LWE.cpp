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

void add_msg(LweSample *result, LweSample *input, const LweParams *lweparams, int mu){
	int n = lweparams->n;
	for(int i = 0; i < n; i++){
		result->a[i] = input->a[i];
	}
	result->b = input->b + modSwitchToTorus32(mu, 2);
}

void generate_sample(LweSample *sample, LweSample *input, const LweParams *lweparams, const TFheGateBootstrappingCloudKeySet *cloud_key, int bit){
	LweSample *xor_res = new_LweSample(lweparams);

	bootsXOR(xor_res, input, input, cloud_key);

	if(bit == 0){
		add_msg(sample, xor_res, lweparams, 0);
	}
	else{
		add_msg(sample, xor_res, lweparams, 1);
	}
}

int main(){
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dist(0, 1);

	auto myparams = initialize_gate_bootstrapping_params();
	auto params = myparams->bootparams;
	const LweParams *lweparams = params->in_out_params;

	TFheGateBootstrappingSecretKeySet *key = new_random_gate_bootstrapping_secret_keyset(params);
	auto cloud_key = &key->cloud;

	/* choose a random bit */
	int32_t sel = dist(gen);

	/* Encryption : Encrypt bit 0 if sel = 0, else encrypt bit 1 */
	LweSample *ciphertext = new_LweSample(lweparams);
	if(sel == 0){		
		bootsSymEncrypt(ciphertext, 0, key);
	}
	else{
		bootsSymEncrypt(ciphertext, 1, key);
	}
	std::cout << "passing encryption of: " << sel << "\n\n";

	/* Generate fresh samples from two distributions, first from the encryptions of 0 and second from the encryptions of 1 */
	LweSample *fresh_sample0 = new_LweSample(lweparams);
	LweSample *fresh_sample1 = new_LweSample(lweparams);

	/* verify correctness */
	std::cout << "Generating sample of encryption of 0\n";
	generate_sample(fresh_sample0, ciphertext, lweparams, cloud_key, 0);
	int dmsg0 = bootsSymDecrypt(fresh_sample0, key);
	std::cout << "expected: 0, generated: " << dmsg0 << "\n\n";

	std::cout << "Generating sample of encryption of 1\n";
	generate_sample(fresh_sample1, ciphertext, lweparams, cloud_key, 1);
	int dmsg1 = bootsSymDecrypt(fresh_sample1, key);
	std::cout << "expected: 1, generated: " << dmsg1 << "\n\n";
	
	return 0;
}