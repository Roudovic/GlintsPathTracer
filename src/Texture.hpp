//
//  Texture.hpp
//  PathTracer
//
//  Created by Ludovic Theobald on 05/11/2019.
//  Copyright © 2019 Ludovic Theobald. All rights reserved.
//

#ifndef Texture_hpp
#define Texture_hpp

#include <stdio.h>
#include <string>
#include "RayVec.hpp"
#include "SrgbTransform.hpp"

class Texture{
public:

    Texture(const std::string & filename, int type, float _scale_u, float _scale_v){
        data = loadTextureFromFile(filename,type);
        
        scale_u = _scale_u;
        scale_v = _scale_v;
    }
    float * loadTextureFromFile(const std::string & filename, int type);
    template<typename T> T queryTexture(float u, float v);
    int get_width() const;
    int get_height() const;
    float get_scale_u() const;
    float get_scale_v() const;
    ~Texture();
    
    
private:
    float * data;
    int width, height, numComponents;
    float scale_u, scale_v;
    
};

#endif /* Texture_hpp */
