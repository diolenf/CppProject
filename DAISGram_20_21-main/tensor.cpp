#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "dais_exc.h"
#include "tensor.h"

#define PI 3.141592654
#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */

using namespace std;


/**
 * Random Initialization
 * 
 * Perform a random initialization of the tensor
 * 
 * @param mean The mean
 * @param std  Standard deviation
 */
void Tensor::init_random(float mean, float std){
    if(data){

        std::default_random_engine generator;
        std::normal_distribution<float> distribution(mean,std);

        for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                for(int k=0;k<d;k++){
                    this->operator()(i,j,k)= distribution(generator);
                }
            }
        }    

    }else{
        throw(tensor_not_initialized());
    }
}

Tensor::Tensor(int r, int c, int d, float v) {
	this->r = r;
	this->c = c;
	this->d = d;
	data = new float**[r];
	for (int i = 0; i < r; i++) {
		data[i] = new float* [c];

		for (int j = 0; j < c; j++) {
			data[i][j] = new float[d];
			for (int k = 0; k < d; k++) {
				data[i][j][k] = v;
			}
		}
	}
}

Tensor::~Tensor(){
	for(int i = 0; i < r; i++){
		for(int j = 0; j < c; j++){
			delete[] data[i][j];
		}
		delete[] data[i];
	}
	delete[] data;
}
