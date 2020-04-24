//
//  shape.hpp
//  PathTracer
//
//  Created by Ludovic Theobald on 13/04/2020.
//  Copyright © 2020 Ludovic Theobald. All rights reserved.
//

#ifndef shape_hpp
#define shape_hpp

#include <stdio.h>
#include "brdf.hpp"
#include "Intersection.hpp"
#include "Texture.hpp"


enum Refl_t { DIFFUSE, MIRROR, GLASS , EMMISSIVE };  // material types, used in radiance()

class Shape{
public:
    Shape( Vec e, Vec c, Refl_t refl,
          bool albedo_textured = false,
          bool normal_textured = false,
          std::shared_ptr<Texture> albedo_tex = nullptr,
          std::shared_ptr<Texture> normal_tex = nullptr,
          std::shared_ptr<BRDF> brdf = std::make_shared<Lambertian>(1)):
    e_(e), c_(c), refl_(refl),
    albedo_textured_(albedo_textured),
    normal_textured_(normal_textured),
    albedo_tex_(albedo_tex),
    normal_tex_(normal_tex),
    brdf_(brdf){
        rayDiff_ = (brdf->getName() == "PatchGGX");
    }
    virtual double intersect(const RayDifferentials &r, Intersection &isect, bool rayDiff) const  = 0;
    virtual Vec randomSample() const = 0;
    virtual Vec randomSample(Vec p_from, double &pdf_li) const = 0;
    
    virtual Vec getEmission() const { return e_;}
    virtual Vec getColor() const { return c_;}
    virtual Refl_t getReflectionType() const { return refl_;}
    virtual bool isAlbedoTextured() const { return albedo_textured_;}
    virtual bool isNormalTextured() const { return normal_textured_;}
    virtual std::shared_ptr<Texture> getAlbedoTex() const{ return albedo_tex_;}
    virtual std::shared_ptr<Texture> getNormalTex() const{ return normal_tex_;}
    virtual std::shared_ptr<BRDF> getBrdf() const { return brdf_;}
    virtual bool needRayDiff() const { return rayDiff_;}


    
protected :
    Vec e_, c_;
    Refl_t refl_;
    bool albedo_textured_, normal_textured_;
    std::shared_ptr<Texture> albedo_tex_, normal_tex_;
    std::shared_ptr<BRDF> brdf_;
    bool rayDiff_;
    
};


#endif /* shape_hpp */
