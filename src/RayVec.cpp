//
//  RayVec.cpp
//  PathTracer
//
//  Created by Ludovic Theobald on 06/11/2019.
//  Copyright © 2019 Ludovic Theobald. All rights reserved.
//

#include <stdio.h>
#include "RayVec.hpp"

Vec operator*(double b , Vec const & o) { return Vec(o.x*b,o.y*b,o.z*b); }


void generateRandomPointOnSphere(double & theta , double & phi) {
    double x = (double)(rand()) / RAND_MAX;
    double y = (double)(rand()) / RAND_MAX;
    theta = x * 2.0 * M_PI;
    phi = acos( std::min<double>(1.0 , std::max<double>(-1.0 , 2.0*y - 1.0 ) ) );
}
Vec randomSampleOnSphere() {
    double theta , phi;
    generateRandomPointOnSphere(theta , phi);
    return Vec( cos(theta)*cos(phi) , sin(theta)*cos(phi) , sin(phi) );
}
Vec randomSampleOnHemisphere( Vec const & upDirection ) {
    Vec r = randomSampleOnSphere();
    if( r.dot(upDirection) > 0.0 ) return r;
    return -1.0 * r;
}
Vec cosineSampleHemisphere(Vec const &upDirection){
    double x = (double)(rand()) / RAND_MAX;
    double y = (double)(rand()) / RAND_MAX;
    double phi = asin(sqrt(x));
    double theta = 2*M_PI*y;
    Vec n = upDirection;
    n.normalize();
    Vec ex = Vec(0.5,0.5,0);
    Vec t = n.cross(ex);
    Vec b = n.cross(t);
    return  cos(theta)*cos(phi)*b + sin(theta)*cos(phi)*t + sin(phi) * n ;
    
}


