//
//  triangle.cpp
//  PathTracer
//
//  Created by Ludovic Theobald on 13/04/2020.
//  Copyright © 2020 Ludovic Theobald. All rights reserved.
//

#include "triangle.hpp"



double Triangle::intersect(const RayDifferentials &r, Intersection &isect, bool rayDiff) const{
    // Translate vertices based on ray origin
    Vec p0t = p0_ - (r.o);
    Vec p1t = p1_ - (r.o);
    Vec p2t = p2_ - (r.o);
    
    
    
    // Apply shear transformation to translated vertex positions
    Vec d = r.d;
    double Sx = -d.x / d.z;
    double Sy = -d.y / d.z;
    double Sz = 1.f / d.z;
    p0t.x += Sx * p0t.z;
    p0t.y += Sy * p0t.z;
    p1t.x += Sx * p1t.z;
    p1t.y += Sy * p1t.z;
    p2t.x += Sx * p2t.z;
    p2t.y += Sy * p2t.z;
    // Compute edge function coefficients _e0_, _e1_, and _e2_
    double e0 = p1t.x * p2t.y - p1t.y * p2t.x;
    double e1 = p2t.x * p0t.y - p2t.y * p0t.x;
    double e2 = p0t.x * p1t.y - p0t.y * p1t.x;
    
    // Perform triangle edge and determinant tests
    if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
        return 0.0;
    double det = e0 + e1 + e2;
    if (det == 0) return 0;
    
    p0t.z *= Sz;
    p1t.z *= Sz;
    p2t.z *= Sz;
    double tScaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
    
    
    
    // Compute barycentric coordinates and $t$ value for triangle intersection
    double invDet = 1 / det;
    double b0 = e0 * invDet;
    double b1 = e1 * invDet;
    double b2 = e2 * invDet;
    double t = tScaled * invDet;
    if(t< 0 ){
        return 0;
    }
    
    Vec edge1 = p1_ - p0_;
    Vec edge2 = p2_ - p0_;
    Vec n =  ((edge1).cross(edge2)).normalize();
    
    //Backface culling
    if(r.d.dot(n) > 0){
//        n = -1.0*n;
        return INFINITY;
    }
//    float orientation = r.d.dot(n) > 0 ? 1 : -1.0;
//    n = orientation * n;
    Vec p_isect = b0*p0_ + b1*p1_ + b2*p2_;
    
    Vec2 uv =  uv0_ * b0 +  uv1_* b1 +  uv2_ * b2;
    
    Vec dpdu, dpdv;
    Vec2 duv02 = uv0_ - uv2_, duv12 = uv1_ - uv2_;
    Vec dp02 = p0_ - p2_, dp12 = p1_ - p2_;
    double determinant = duv02.x * duv12.y - duv02.y * duv12.x;
    if (determinant == 0) {
    } else {
        float invdet = 1 / determinant;
        dpdu = ( duv12.y * dp02 - duv02.y * dp12) * invdet;
        dpdv = (-duv12.x * dp02 + duv02.x * dp12) * invdet;
    }
    
    
    isect = Intersection(p_isect,n, uv.x, uv.y, dpdu,dpdv, 0);
    if(rayDiff){
        isect.computeFootprint(r, 1./(32));
    }
//    float dist = (p_isect - r.o).norm();

    return t;
}
Vec Triangle::randomSample() const {
    double a, b,  b0, b1, b2;
    a = (double)rand()/RAND_MAX;
    b = (double)rand()/RAND_MAX;
    if(a>b){
        std::swap(a,b);
    }
    b0 = a;
    b1 = b - a;
    b2 = 1 -b;
    
    return  b0*p0_ + b1*p1_ + b2*p2_;
}
Vec Triangle::randomSample(Vec p_from, double &pdf_li) const {
    Vec p_samp = randomSample();
    Vec edge1 = p1_ - p0_;
    Vec edge2 = p2_ - p0_;
    double area = ((edge1).cross(edge2)).norm();
    double dist = (p_samp - p_from).norm();
    pdf_li = 4*M_PI / ( area/(dist*dist) );
    
    return (p_samp - p_from).normalize();
    
    //TODO : Compute pdf_li
}
