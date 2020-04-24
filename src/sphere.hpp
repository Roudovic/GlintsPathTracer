//
//  sphere.hpp
//  PathTracer
//
//  Created by Ludovic Theobald on 13/04/2020.
//  Copyright © 2020 Ludovic Theobald. All rights reserved.
//

#ifndef sphere_h
#define sphere_h

#include "shape.hpp"

class Sphere : public Shape {
private:
    double radius;      // radius
    Vec p_;
public:
    Sphere( double rad, Vec p, Vec e, Vec c, Refl_t refl,
           bool albedo_textured = false,
           bool normal_textured = false,
           std::shared_ptr<Texture> albedo_tex = nullptr,
           std::shared_ptr<Texture> normal_tex = nullptr,
           std::shared_ptr<BRDF> brdf = std::make_shared<Lambertian>(1)):
    radius(rad),p_(p), Shape(e,c, refl, albedo_textured, normal_textured, albedo_tex, normal_tex, brdf){}
    virtual double handle_intersect(const double lambda, const RayDifferentials &r, Intersection &isect, bool rayDiff) const ;
    virtual double intersect(const RayDifferentials &r, Intersection &isect, bool rayDiff) const ;
    virtual Vec randomSample() const ;
    virtual Vec randomSample(Vec p_from, double &pdf_li) const ;
    
    virtual Vec getPosition() const { return p_;}

};

#endif /* sphere_h */
