//
//  mesh.hpp
//  PathTracer
//
//  Created by Ludovic Theobald on 18/04/2020.
//  Copyright © 2020 Ludovic Theobald. All rights reserved.
//

#ifndef mesh_hpp
#define mesh_hpp

#include <stdio.h>
#include "shape.hpp"
#include "triangle.hpp"
#include "sphere.hpp"
class Mesh : public Shape{
public:
    Mesh(const std::string& filename, Vec translation, Mat3 rotation, double scale, Vec e, Vec c, Refl_t refl,
         bool albedo_textured = false,
         bool normal_textured = false,
         std::shared_ptr<Texture> albedo_tex = nullptr,
         std::shared_ptr<Texture> normal_tex = nullptr,
         std::shared_ptr<BRDF> brdf = std::make_shared<Lambertian>(1)):
    Shape(e,c, refl, albedo_textured, normal_textured, albedo_tex, normal_tex, brdf){
        loadOFF(filename, translation, rotation, scale);
    }
    void loadOFF(const std::string& filename, Vec translation, Mat3 rotation, double scale);
    void applyTransform(Vec translation, Mat3 rotation, double scale);
    virtual double intersect(const RayDifferentials &r, Intersection &isect, bool rayDiff) const ;
    virtual Vec randomSample() const ;
    virtual Vec randomSample(Vec p_from, double &pdf_li) const ;


    
private:
    std::vector<Triangle> triangles_;
    std::unique_ptr<Sphere> boundingSphere_;
};
#endif /* mesh_hpp */
