#include <iostream>
#include <string>

#include "dais_exc.h"
#include "tensor.h"
#include "libbmp.h"
#include "DAISGram.h"

using namespace std;

DAISGram::DAISGram() {}

DAISGram::~DAISGram() {}

/**
 * Load a bitmap from file
 *
 * @param filename String containing the path of the file
 */
void DAISGram::load_image(string filename){
    BmpImg img = BmpImg();

    img.read(filename.c_str());

    const int h = img.get_height();
    const int w = img.get_width();

    data = Tensor(h, w, 3, 0.0);

    for(int i=0;i<img.get_height();i++){
        for(int j=0;j<img.get_width();j++){ 
            data(i,j,0) = (float) img.red_at(j,i);
            data(i,j,1) = (float) img.green_at(j,i);    
            data(i,j,2) = (float) img.blue_at(j,i);   
        }                
    }
}


/**
 * Save a DAISGram object to a bitmap file.
 * 
 * Data is clamped to 0,255 before saving it.
 *
 * @param filename String containing the path where to store the image.
 */
void DAISGram::save_image(string filename){

    data.clamp(0,255);

    BmpImg img = BmpImg(getCols(), getRows());

    img.init(getCols(), getRows());

    for(int i=0;i<getRows();i++){
        for(int j=0;j<getCols();j++){
            img.set_pixel(j,i,(unsigned char) data(i,j,0),(unsigned char) data(i,j,1),(unsigned char) data(i,j,2));                   
        }                
    }

    img.write(filename);

}


/**
 * Generate Random Image
 * 
 * Generate a random image from nois
 * 
 * @param h height of the image
 * @param w width of the image
 * @param d number of channels
 * @return returns a new DAISGram containing the generated image.
 */  
void DAISGram::generate_random(int h, int w, int d){
    data = Tensor(h,w,d,0.0);
    data.init_random(128,50);
    data.rescale(255);
}


DAISGram DAISGram::grayscale(){
    float average=0.0;
    DAISGram output;
    output.data.init(getRows(),getCols(),getDepth());
    for(int i=0;i<getRows();i++){
        for(int j=0;j<getCols();j++){
            for(int k=0;k<getDepth();k++){
                average=average+data(i,j,k);
            }
            average=average/getDepth();
            for(int l=0;l<output.getDepth();l++){
                output.data(i,j,l)=average;
            }
            average=0.0;
        }
    }
    return output;
}

DAISGram DAISGram::brighten(float bright){
    DAISGram output;
    output.data = this->data;
    output.data = output.data + bright;
    output.data.clamp(0,255);
    return output;
}



/*tested sharpen */ 
DAISGram DAISGram::sharpen(){
    DAISGram output;
    Tensor filter;
    filter.read_file("Filters//sharpen.txt");
    output.data=data.convolve(filter);
    output.data.clamp(0,255);
    return output;
}

DAISGram DAISGram::emboss(){
    DAISGram output;
    Tensor filter;
    filter.read_file("Filters//emboss.txt");
    output.data=data.convolve(filter);
    output.data.clamp(0,255);
    return output;
}

DAISGram DAISGram::smooth(int h){
    float c=1.f/((float)h*h);
    Tensor filter(h,h,1,c);
    DAISGram output;
    output.data=data.convolve(filter);
    return output;
}

DAISGram DAISGram::edge(){
    Tensor filter;
    DAISGram output;
    filter.read_file("Filters//edge.txt");
    output=grayscale();
    output.data=output.data.convolve(filter);
    output.data.clamp(0,255);
    return output;
}



Tensor DAISGram::swap_channel(int index1, int index2){
	Tensor res = data;
	for (int i = 0; i < res.rows(); i++) {
		for (int j = 0; j < res.cols(); j++) {
			float aux = res(i, j, index1);
			res(i, j, index1) = res(i, j, index2);
			res(i, j, index2) = aux;
		}
	}
	return res;
}

DAISGram DAISGram::warhol(){
	DAISGram res;
	Tensor right_up = swap_channel(0, 1);       //Red <=> green
	Tensor left_bot = swap_channel(1, 2);      //Blue <=> Green
	Tensor right_bot = swap_channel(0, 2);    //Red <=> Blue
	Tensor concat_up = data.concat(right_up,1);
	Tensor concat_bot = left_bot.concat(right_bot,1);
	Tensor last = concat_up.concat(concat_bot);
	res.data = last;
	return res;
}

int DAISGram::getRows(){
	return data.rows();
}

int DAISGram::getCols() {
	return data.cols();
}

int DAISGram::getDepth() {
	return data.depth();
}


DAISGram DAISGram::blend(const DAISGram &rhs, float alpha) { 
	if(alpha <0 || alpha>1) 
		throw unknown_exception(); 
	DAISGram res; 
	res.data = (data *alpha) + (rhs.data * (1-alpha)); 
	return res;
}


DAISGram DAISGram::greenscreen(DAISGram& bkg, int rgb[], float threshold[]){
	if (getRows() != bkg.getRows() || getCols() != bkg.getCols() || getDepth() != bkg.getDepth())
		throw dimension_mismatch();
	DAISGram res;
	res.data = this->data;
	bool red, green, blue;
	red = green = blue = false;
	for (int i = 0; i < getRows(); i++) {
		for (int j = 0; j < getCols(); j++) {
			if (res.data(i, j, 0) >= rgb[0] - threshold[0] && res.data(i, j, 0) <= rgb[0] + threshold[0])
				red = true;
			if (res.data(i, j, 1) >= rgb[1] - threshold[1] && res.data(i, j, 1) <= rgb[1] + threshold[1])
				green = true;
			if (res.data(i, j, 2) >= rgb[2] - threshold[2] && res.data(i, j, 2) <= rgb[2] + threshold[2])
				blue = true;
			if (red && green && blue) {
				for (int k = 0; k < getDepth(); k++) {
					res.data(i, j, k) = bkg.data(i, j, k);
				}
			}
			red = blue = green = false;
		}
	}
	return res;
}

DAISGram DAISGram::equalize(){
	throw method_not_implemented();
}