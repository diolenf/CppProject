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

Tensor::Tensor(const Tensor& that) {
	r = that.r;
	c = that.c;
	d = that.d;
	data = new float** [r];
	for (int i = 0; i < r; i++) {
		data[i] = new float* [c];

		for (int j = 0; j < c; j++) {
			data[i][j] = new float[d];
		}
	}
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			for (int k = 0; k < d; k++) {
				data[i][j][k] = that.data[i][j][k];
			}
		}
	}

}

bool Tensor::operator==(const Tensor& rhs) const {
	if (c != rhs.c || r != rhs.r || d != rhs.d)
		throw dimension_mismatch();
	else if(!data)
		throw tensor_not_initialized();
	else {
		bool rs = true;
		int i = 0;
		while (i < r && rs) {
			int j = 0;
			while (j < c && rs) {
				int k = 0;
				while (k < d && rs) {
					if (std::fabs(data[i][j][k] - rhs.data[i][j][k]) >= EPSILON) {
						rs = false;
					}
					k++;
				}
				j++;
			}
			i++;
		}
		return rs;
	}
}

int Tensor::rows() const{
	return r;
}

int Tensor::cols() const{
	return c;
}

int Tensor::depth() const{
	return d;
}

float Tensor::getMin(int k) const{
	if(!data)
		throw tensor_not_initialized();
	
	float min = this->operator()(0, 0, k);
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			if (this->operator()(i, j, k) < min)
				min = this->operator()(i, j, k);
		}
	}
	return min;
}

float Tensor::getMax(int k) const{
	if(!data)
		throw tensor_not_initialized();
	
	float max = this->operator()(0, 0, k);
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			if (this->operator()(i, j, k) > max)
				max = this->operator()(i, j, k);
		}
	}
	return max;
}

void Tensor::showSize() const{
	std::cout << r << "x" << c << "x" << d;
}

std::ostream& operator<<(std::ostream& stream, const Tensor& obj){
	if(!obj.data)
		stream << "[]";
	else {
		for (int i = 0; i < obj.r; i++) {
			for (int j = 0; j < obj.c; j++) {
				stream << "[";
				for (int k = 0; k < obj.d; k++) {
					stream << obj(i, j, k) << ",";
				}
				stream << "\b] ";
			}
			stream << std::endl;
		}
	}
	return stream;
}

Tensor Tensor::operator-(const Tensor& rhs)const{
  Tensor a;
 a.init(r,c,d,0.0);
  for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			for (int k = 0; k < d; k++) {
				a.data(i,j,k) = data[i][j][k]-rhs.data(i,j,k);
			}
		}
	}
	return a;
}

Tensor Tensor::operator+(const Tensor& rhs)const{
   Tensor a;
 a.init(r,c,d,0.0);
  for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			for (int k = 0; k < d; k++) {
				a.data(i,j,k) = data[i][j][k]+rhs.data(i,j,k);
			}
		}
	}
	return a;
}
Tensor Tensor::operator*(const Tensor& rhs)const{
   Tensor a;
 a.init(r,c,d,0.0);
  for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			for (int k = 0; k < d; k++) {
				a.data(i,j,k) = data[i][j][k]*rhs.data(i,j,k);
			}
		}
	}
	return a;
}
Tensor Tensor::operator/(const Tensor& rhs)const{
  Tensor a;
 a.init(r,c,d);
  for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			for (int k = 0; k < d; k++) {
				a.data(i,j,k) = data[i][j][k]/rhs.data(i,j,k);
			}
		}
	}
	return a;
}

Tensor Tensor::operator-(const float& rhs)const{
  Tensor a;
  a.init(r,c,d);
  for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			for (int k = 0; k < d; k++) {
				a.data(i,j,k) = data[i][j][k]-rhs;
			}
		}
	}
	return a;
}

Tensor Tensor::operator+(const float& rhs)const{
  Tensor a;
  a.init(r,c,d);
  for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			for (int k = 0; k < d; k++) {
				a.data(i,j,k) = data[i][j][k]+rhs;
			}
		}
	}
	return a;
}

Tensor Tensor::operator*(const float& rhs)const{
 Tensor a;
  a.init(r,c,d);
  for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			for (int k = 0; k < d; k++) {
				a.data(i,j,k) = data[i][j][k]*rhs;
			}
		}
	}
	return a;
}

Tensor Tensor::operator/(const float& rhs)const{
 Tensor a;
  a.init(r,c,d);
  for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			for (int k = 0; k < d; k++) {
				a.data(i,j,k) = data[i][j][k]/rhs;
			}
		}
	}
	return a;
}

void Tensor::clamp(float low,float high){
   for(int i=0;i<r;i++){
     for(int j=0;j<c;j++){
	   for(int k=0;k<d;k++){
	     if(data[i][j][k])<low){
		 data[i][j][k]=low;
		  }
		 if(data[i][j][k]>high){
		  data[i][j][k]=high;
		  }
		}
	  }
	}
}




Tensor Tensor::padding(int pad_h, int pad_w)const {
	
	if(!data)
		throw tensor_not_initialized();
	
	Tensor aux( r + 2 * pad_h, c + 2 * pad_w, d);  
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < c; j++) {
			for (int k = 0; k < d; k++) {
				aux(i + pad_h, j + pad_w, k) = data[i][j][k];
			}
		}
	}
	return aux;
}

Tensor Tensor::subset(unsigned int row_start, unsigned int row_end, unsigned int col_start,
			 unsigned int col_end, unsigned int depth_start, unsigned int depth_end) const
{
	
	if(!data)
		throw tensor_not_initialized();
	
	int new_r = row_end - row_start;
	int new_c = col_end - col_start;
	int new_d = depth_end - depth_start;
	
	if(new_r <= 0 || new_c <= 0 || new_d <= 0)
		throw unknown_exception();

	Tensor newt(new_r, new_c, new_d);
	for (int i = row_start; i < row_end; i++) {
		for (int j = col_start; j < col_end; j++) {
			for (int k = depth_start; k < depth_end; k++) {
				newt(i - row_start, j - col_start, k - depth_start) = this->operator()(i,j,k);
			}
		}
	}
	return newt;
}


void Tensor::rescale(float new_max=1.0){
    for(int i=0;i<r;i++){
        for(int j=0;j<c;j++){
            for(int k=0;k<d;k++){
                if(getMax(k)!=getMin(k)){
                    data[i][j][k]=((data[i][j][k]-getMin(k))/(getMax(k)-getMin(k)))*new_max;
                }
                else{
                    data[i][j][k]=new_max;
                }
            }
        }
    }
}

void Tensor::read_file(string filename){
	std::ifstream file(filename);
	if (!file.is_open())
		throw unable_to_read_file();
	file >> r >> c >> d;
	init(r, c, d);
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < r; j++) {
			for (int k = 0; k < c; k++) {
				float temp;
				file >> temp;
				this->operator()(j, k, i) = temp;
			}
		}
	}
	file.close();
}

void Tensor::write_file(string filename) {
	std::ofstream file(filename);
	file << r << "\n";
	file << c << "\n";
	file << d << "\n";
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < r; j++) {
			for (int k = 0; k < c; k++) {
				file << this->operator()(j, k, i) <<"\n";
			}
		}
	}
	file.close();
}

Tensor Tensor::convolve(const Tensor& f)const {
	if (!data)
		throw tensor_not_initialized();
	if (f.r % 2 != 1)
		throw filter_odd_dimensions();
	if (f.r != f.c)
		throw unknown_exception();

	int pad = (f.r - 1) / 2;
	Tensor res(r, c, d, 0.0);
	Tensor padded = padding(pad, pad);

	for (int i = 0; i <= padded.r - f.r; i++) {
		for (int j = 0; j <= padded.c - f.c; j++) {
			for (int k = 0; k < d; k++) {
				Tensor aux = padded.subset(i, i + f.r, j, j + f.c, k, k + 1); // 
				float sum = conv_sum(aux, f);
				res(i, j, k) = sum;
			}
		}
	}
	return res;
}


