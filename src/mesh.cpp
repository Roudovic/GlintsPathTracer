//
//  mesh.cpp
//  PathTracer
//
//  Created by Ludovic Theobald on 18/04/2020.
//  Copyright © 2020 Ludovic Theobald. All rights reserved.
//

#include "mesh.hpp"
#include <iostream>
#include <fstream>
#include <exception>
#include <ios>
#include <vector>
#include <array>

void Mesh::loadOFF (const std::string & filename, Vec translation, Mat3 rotation, double scale) {
    std::cout << " > Start loading mesh <" << filename << ">" << std::endl;
    std::ifstream in (filename.c_str ());
    if (!in)
        throw std::ios_base::failure ("[Mesh Loader][loadOFF] Cannot open " + filename);
    std::string offString;
    unsigned int sizeV, sizeT, tmp;
    in >> offString >> sizeV >> sizeT >> tmp;
//    auto & P = meshPtr->vertexPositions ();
//    auto & T = meshPtr->triangleIndices ();
    std::vector<std::array<double, 3>> P(sizeV);
    std::vector<std::array<int, 3>> T(sizeT);

    size_t tracker = (sizeV + sizeT)/20;
    std::cout << " > [" << std::flush;
    for (unsigned int i = 0; i < sizeV; i++) {
        if (i % tracker == 0)
            std::cout << "-" << std::flush;
        in >> P[i][0] >> P[i][1] >> P[i][2];
    }
    int s;
    for (unsigned int i = 0; i < sizeT; i++) {
        if ((sizeV + i) % tracker == 0)
            std::cout << "-" << std::flush;
        in >> s;
        for (unsigned int j = 0; j < 3; j++)
            in >> T[i][j];
    }
    std::cout << "]" << std::endl;
    in.close ();
    
    for(int i = 0; i < sizeT; i++){
        Vec p0, p1, p2;
//        p0 =  Vec(P[T[i][0]][0],P[T[i][0]][1],P[T[i][0]][2]) + translation;
//        p1 =  Vec(P[T[i][1]][0],P[T[i][1]][1],P[T[i][1]][2]) + translation;
//        p2 =  Vec(P[T[i][2]][0],P[T[i][2]][1],P[T[i][2]][2]) + translation;
        
        
        p0 =  scale *(rotation.mult(Vec(P[T[i][0]][0],P[T[i][0]][1],P[T[i][0]][2]))) + translation;
        p1 =  scale *(rotation.mult(Vec(P[T[i][1]][0],P[T[i][1]][1],P[T[i][1]][2]))) + translation;
        p2 =  scale *(rotation.mult(Vec(P[T[i][2]][0],P[T[i][2]][1],P[T[i][2]][2]))) + translation;
        
        triangles_.emplace_back(p0,p1,p2,
                                e_, c_, refl_,
                                albedo_textured_,
                                normal_textured_,
                                albedo_tex_,
                                normal_tex_,
                                brdf_);
        
    }
    double max_radius = 0;
    for(int i = 0; i < sizeV; i++){
        Vec p0 = Vec(P[i][0], P[i][1], P[i][2]);
        double radius = p0.norm();
        if(radius > max_radius){
            max_radius = radius;
        }
    }
    boundingSphere_ = make_unique<Sphere>(scale*max_radius, translation, e_, c_, refl_,
                             albedo_textured_,
                             normal_textured_,
                             albedo_tex_,
                             normal_tex_,
                             brdf_);
    //    meshPtr->recomputePerVertexNormals ();
    std::cout << " > Mesh <" << filename << "> loaded" <<  std::endl;
}

double Mesh::intersect(const RayDifferentials &r, Intersection &isect, bool rayDiff) const{
    float dist = INFINITY;
    if(boundingSphere_->intersect(r, isect, false) > 0){
        for (auto& triangle : triangles_){
            float cur_dist = triangle.intersect(r, isect, false);
            if( cur_dist > 0 && cur_dist < dist){
                dist = cur_dist;
            }
        }
        if(rayDiff){
            isect.computeFootprint(r, 1./16);
        }
        return dist;
    }
    else{
        return 0;
    }
}
Vec Mesh::randomSample() const {
    int rand_id = rand() % triangles_.size();
    return triangles_[rand_id].randomSample();
}
Vec Mesh::randomSample(Vec p_from, double &pdf_li) const {
    Vec p_samp = randomSample();
    //TODO : Figure out the pdf
    pdf_li = 4*M_PI / 100;
    
    return (p_samp - p_from).normalize();
}
