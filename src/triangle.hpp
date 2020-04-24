//
//  triangle.hpp
//  PathTracer
//
//  Created by Ludovic Theobald on 13/04/2020.
//  Copyright © 2020 Ludovic Theobald. All rights reserved.
//

#ifndef triangle_hpp
#define triangle_hpp

#include <stdio.h>
#include "shape.hpp"


class Triangle : public Shape {
private:
    Vec p0_, p1_, p2_;
    Vec2 uv0_, uv1_, uv2_;
public:
    Triangle(Vec p0, Vec p1, Vec p2, Vec e, Vec c, Refl_t refl,
           bool albedo_textured = false,
           bool normal_textured = false,
           std::shared_ptr<Texture> albedo_tex = nullptr,
           std::shared_ptr<Texture> normal_tex = nullptr,
           std::shared_ptr<BRDF> brdf = std::make_shared<Lambertian>(1),
             Vec2 uv0 =  Vec2(0,0), Vec2 uv1 = Vec2(1,0), Vec2 uv2 =Vec2(1,1)):
    p0_(p0), p1_(p1),p2_(p2), Shape(e,c, refl, albedo_textured, normal_textured, albedo_tex, normal_tex, brdf), uv0_(uv0), uv1_(uv1), uv2_(uv2){    
    }
    virtual double intersect(const RayDifferentials &r, Intersection &isect, bool rayDiff) const ;
    virtual Vec randomSample() const ;
    virtual Vec randomSample(Vec p_from, double &pdf_li) const ;
    
    virtual Vec getP0() const { return p0_;}
    virtual Vec getP1() const { return p1_;}
    virtual Vec getP2() const { return p2_;}


};
#endif /* triangle_hpp */
