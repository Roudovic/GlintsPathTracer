//
//  Intersection.hpp
//  PathTracer
//
//  Created by Ludovic Theobald on 20/11/2019.
//  Copyright © 2019 Ludovic Theobald. All rights reserved.
//

#ifndef Intersection_hpp
#define Intersection_hpp

#include <stdio.h>
#include "RayVec.hpp"

class Intersection{
public:
    Intersection() = default;

    Intersection(Vec p_isect, Vec n,float u,float v, Vec dpdu, Vec dpdv, Mat2 footprint, int depth=0):
    p_isect(p_isect),
    ng(n),
    n(n),
    u(u), v(v),
    dpdu(dpdu),dpdv(dpdv),
    footprint(footprint),
    depth(depth){}
    float get_u() const{
        return u;
    }
    float get_v() const{
        return v;
    }
    Vec get_n() const{
        return n;
    }
    void set_n(const Vec& n_){
        n = n_;
    }

    Vec get_ng() const{
        return ng;
    }
    Vec get_dpdu() const{
        return dpdu;
    }
    Vec get_dpdv() const{
        return dpdv;
    }
    Vec get_dpdx() const{
        return dpdx;
    }
    Vec get_dpdy() const{
        return dpdy;
    }
    Mat2 get_fp() const{
        return footprint;
    }
    int get_depth() const{
        return depth;
    }
    void set_depth(int d){
        depth=d;
    }
    void computeFootprint(const RayDifferentials& r, float scale_fp = 1./(1024*16));
private:
    Vec p_isect, ng, n;
    float u,v;
    Vec dpdu, dpdv, dpdx, dpdy;
    Mat2 footprint;
    int depth;
};
#endif /* Intersection_hpp */
