//
//  position_normal_bvh.hpp
//  PathTracer
//
//  Created by Ludovic Theobald on 10/01/2020.
//  Copyright © 2020 Ludovic Theobald. All rights reserved.
//

#ifndef position_normal_bvh_hpp
#define position_normal_bvh_hpp

#include <stdio.h>
#include <stdlib.h>
#include "Texture.hpp"
#include "RayVec.hpp"
#include "position_normal_bvh.hpp"

class BoundingBox{
public:
    BoundingBox(Vec2 v):min_(v), max_(v){}
    BoundingBox():min_(0), max_(0){}
private:
    
    Vec2 min_, max_;
    
};


class PositionNormalBvh{
public:
    PositionNormalBvh(double intrinsic_roughness, BoundingBox bbox): intrinsic_roughness_(intrinsic_roughness), bbox_(bbox){};
    std::vector<Vec2> intersectBvh(double min_u, double min_v, double max_u, double max_v, Vec2 s);
    
private:
    double intrinsic_roughness_;
    BoundingBox bbox_;
    PositionNormalBvh *left, *right;
};




#endif /* position_normal_bvh_hpp */

