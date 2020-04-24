//
//  Intersection.cpp
//  PathTracer
//
//  Created by Ludovic Theobald on 20/11/2019.
//  Copyright © 2019 Ludovic Theobald. All rights reserved.
//

#include "Intersection.hpp"

void Intersection::computeFootprint(const RayDifferentials& ray, float scale_fp){
    
    if (ray.hasDifferentials){
        float d = -n.dot(p_isect);
        float tx = (-n.dot(ray.ox) - d)/(n.dot(ray.dir_x));
        Vec px = ray.ox + tx * ray.dir_x;
        float ty = (-n.dot(ray.oy) - d)/(n.dot(ray.dir_y));
        Vec py = ray.oy + ty * ray.dir_y;
        dpdx = px - p_isect;
        dpdy = py - p_isect;
        
        Mat2 A;
        Vec2 bx, by;
        if (std::abs(n.x) > std::abs(n.y) && std::abs(n.x) > std::abs(n.z)) {
            A = Mat2(dpdu.y, dpdu.z, dpdv.y, dpdv.z);
            bx = Vec2(dpdx.y, dpdx.z);
            by = Vec2(dpdy.y, dpdy.z);
        } else if (std::abs(n.y) > std::abs(n.z)) {
            A = Mat2(dpdu.x, dpdu.z, dpdv.x, dpdv.z);
            bx = Vec2(dpdx.x, dpdx.z);
            by = Vec2(dpdy.x, dpdy.z);
            
        } else {
            A = Mat2(dpdu.x, dpdu.y, dpdv.x, dpdv.y);
            bx = Vec2(dpdx.x, dpdx.y);
            by = Vec2(dpdy.x, dpdy.y);
        }
//        A = Mat2(dpdu.y, dpdu.z, dpdv.y, dpdv.z);
//        bx = Vec2(dpdx.y, dpdx.z);
//        by = Vec2(dpdy.y, dpdy.z);


        
        Mat2 A_inv = A.inv();
        Vec2 ddx = A_inv.mult(bx);  // (dudx, dvdx)
        Vec2 ddy = A_inv.mult(by);  // (dudy, dvdy)
        footprint = Mat2(ddx.x, ddy.x, ddx.y, ddy.y) * scale_fp;
       
    }
}
