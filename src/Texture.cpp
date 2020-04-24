//
//  Texture.cpp
//  PathTracer
//
//  Created by Ludovic Theobald on 05/11/2019.
//  Copyright © 2019 Ludovic Theobald. All rights reserved.
//

#include "Texture.hpp"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <math.h>

float * Texture::loadTextureFromFile(const std::string & filename, int type){
// Loading the image in CPU memory using stbd_image
float * data = stbi_loadf (filename.c_str (), &width,
                                  &height,
                                  &numComponents, // 1 for a 8 bit greyscale image, 3 for 24bits RGB image
                                  0);
    
    //stbi_image_free(data);
    std::cout << height << std::endl;
    std::cout << width << std::endl;
    return data;
}


template<> Vec Texture::queryTexture<Vec>(float u, float v){
    int i = (int)floor(u*width*scale_u) % width;
    int j = (int)floor(v*height*scale_v) % height;
    
    float& r = data[numComponents*(width*i+j)];
    float& g = data[numComponents*(width*i+j)+1];
    float& b = data[numComponents*(width*i+j)+2];
    
    return Vec(SrgbTransform::srgbToLinear(r),SrgbTransform::srgbToLinear(g),SrgbTransform::srgbToLinear(b));
}

int Texture::get_width() const{
    return width;
}
int Texture::get_height() const{
    return height;
}
float Texture::get_scale_u() const{
    return scale_u;
}
float Texture::get_scale_v() const{
    return scale_v;
}
Texture::~Texture(){
    stbi_image_free(data);
}
