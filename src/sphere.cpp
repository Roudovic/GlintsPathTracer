
//
//  sphere.cpp
//  PathTracer
//
//  Created by Ludovic Theobald on 13/04/2020.
//  Copyright © 2020 Ludovic Theobald. All rights reserved.
//

#include "sphere.hpp"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include "RayVec.hpp"


double Sphere::handle_intersect(const double lambda, const RayDifferentials &r, Intersection &isect, bool rayDiff) const{
        Vec p_isect = r.o + lambda*r.d;
        Vec ox = p_isect - p_;
        Vec n = ox;
        n.normalize();
        double phi = atan2f(ox.y, ox.x);
        if (phi < 0) phi += 2 * M_PI;
        double cosPhi = cosf(phi);
        double sinPhi = sinf(phi);
        double cosTheta = ox.z/radius;
        double theta = acosf(cosTheta);
        double sinTheta = sinf(theta);
        
        Vec dpdu, dpdv;
        dpdu = Vec(-radius*sinTheta*sinPhi,
                   radius*sinTheta*cosPhi,
                   0);
        dpdv = Vec(radius*cosTheta*cosPhi,
                   radius*cosTheta*cosPhi,
                   -radius*sinTheta);
//        if(rayDiff){
//            // Compute the necessary values for the pixel footprint projection
//            int w = 1024;
//            int h = 768;
//            // The differentials should be aligned with the camera/pellicule axis for this to work.
//            double dx = 1./(w), dy = 1./(h);
//            double dudx, dudy, dvdx,dvdy = 0.0;
//            double phi_dudx = atan2(ox.y, ox.x + dx);
//            double phi_dudy = atan2(ox.y + dy, ox.x);
//            if(phi_dudx < 0) phi_dudx += 2*M_PI;
//            if(phi_dudy < 0) phi_dudy += 2*M_PI;
//            dudx = (1./dx)*(phi_dudx - phi);
//            dudy = (1./dy)*(phi_dudy - phi);
//            int neg_z = theta > M_PI/2 ? -1 : 1;
//            double zdzdx = sqrt(radius*radius - ox.y*ox.y - (ox.x + dx)*(ox.x + dx)) * neg_z;
//            double zdzdy = sqrt(radius*radius - ox.x*ox.x - (ox.y + dy)*(ox.y + dy)) * neg_z;
//            dvdx = (1./dx)*(acos(zdzdx/radius) - theta);
//            dvdy = (1./dy)*(acos(zdzdy/radius) - theta);
//            if(isnan(dvdx))
//                dvdx = 0;
//            if(isnan(dvdy))
//                dvdy = 0;
//
//            //Hard-coded camera related values
//            //TODO : Refactor the code such that the camera parameters are accessible from the path tracer methods
//            double cx = w*.5135/h;
//            double width_pixel = (140*cx/w); //Width of a pixel in world space
//
    //
    //            // fp is the jacobian of (u,v) = w(x,y) x and y being the world space coordinates, scaled by the width of a pixel divided by 16,
    //            // since we implicitly use a gaussian pixel filter whose width is a 16th of a pixel (cf. paper)
    //
    //            Mat2 fp = Mat2(dudx,dudy,dvdx,dvdy) * width_pixel * (1./16);
    //            isect = Intersection(p_isect,n, phi/(2*M_PI), theta/(M_PI), dpdu, dpdv, fp);
    //        }
    //        else{
    //
    //            isect = Intersection(p_isect,n, phi/(2*M_PI), theta/(M_PI), dpdu,dpdv, 0);
    //        }
    isect = Intersection(p_isect, n,phi/(2*M_PI), theta/(M_PI), dpdu,dpdv, 0);
    if(rayDiff){
        isect.computeFootprint(r, 1./(16));
    }
    return lambda;
}

    double Sphere::intersect(const RayDifferentials &r, Intersection &isect, bool rayDiff) const { // returns distance, 0 if nohit
        Vec oc = r.o - p_;
        double sa = 1.0;
        double sb = 2.0 * oc.dot( r.d );
        double sc = oc.dot(oc) - radius * radius;
        
        double delta = sb*sb - 4.0 * sa * sc;
        if( delta < 0.0 )
            return 0.0; // no solution
        
        double deltaSqrt = sqrt(delta);
        double lambda1 = (-sb - deltaSqrt) / (2.0 * sa);
        double lambda2 = (-sb + deltaSqrt) / (2.0 * sa);
        if( lambda1 < lambda2  &&  lambda1 > 0.0 ){
            return handle_intersect(lambda1, r, isect, rayDiff);
        }
        if( lambda2 >= 0.0 ){
            return handle_intersect(lambda2, r, isect, rayDiff);
        }
        return 0.0;
    }
    
    Vec Sphere::randomSample() const {
        return p_ + radius * randomSampleOnSphere();
    }
    
    // Adapted from PBRT v3
    Vec Sphere::randomSample(Vec p_from, double &pdf_li) const{
        double x = (double)(rand()) / RAND_MAX;
        double y = (double)(rand()) / RAND_MAX;
        Vec2 u = Vec2(x,y);
        double sinThetaMax2 = radius * radius / (p_from-p_).dot(p_from-p_);
        double cosThetaMax = sqrt(std::max((double)0, 1 - sinThetaMax2));
        double cosTheta = (1 - u.x) + u.x * cosThetaMax;
        double sinTheta = sqrt(fmax((double)0, 1 - cosTheta * cosTheta));
        double phi = u.y * 2 * M_PI;
        Vec n = (p_ - p_from).normalize();
        Vec ex = Vec(1,0,0);
        Vec t = n.cross(ex).normalize();
        Vec b = n.cross(t).normalize();
        
        pdf_li = 1 / (2*M_PI * (1 - cosThetaMax));
//        pdf_li = (1 - cosThetaMax);
        return cosTheta*n + sin(phi)*sinTheta*t + cos(phi)*sinTheta*b;
        
    }

