//
//  RayVec.hpp
//  PathTracer
//
//  Created by Ludovic Theobald on 05/11/2019.
//  Copyright © 2019 Ludovic Theobald. All rights reserved.
//

#ifndef RayVec_hpp
#define RayVec_hpp

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>




struct Vec {
    double x, y, z;
    Vec(double x_=0, double y_=0, double z_=0){ x=x_; y=y_; z=z_; }
    Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); }
    Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); }
    Vec operator*(double b) const { return Vec(x*b,y*b,z*b); }
    double& operator[](int i){
        if(i==0) return x;
        else if(i==1) return y;
        else return z;

    }
    double operator[](int i) const{
        if(i==0) return x;
        else if(i==1) return y;
        else return z;

    }
    Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); }
    inline Vec& normalize(){ return *this = *this * (1.0/sqrt(x*x+y*y+z*z)); }
    double norm(){ return sqrt(x*x+y*y+z*z);}
    double dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; }
    Vec cross(Vec&b) const{return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}
    Vec TangentSpaceToWorld(Vec n, Vec t, Vec b){
        t.normalize();
        b.normalize();
        //        b = n.cross(t).normalize();
        //        Vec n_ = Vec(this->dot(t),this->dot(b),this->dot(n) ).normalize();
        //        return n_;
        return (n * z + t * x + b * y).normalize();
        //        return t.cross(b) + n.cross(b)*x + t.cross(n)*y;
    }
    Vec toTangentSpace(Vec n, Vec t, Vec b){
        t.normalize();
        b.normalize();
        Vec n_ = Vec(this->dot(t),this->dot(b),this->dot(n) ).normalize();
        return n_;
    }
};

Vec operator*(double b , Vec const & o);


void generateRandomPointOnSphere(double & theta , double & phi);

Vec randomSampleOnSphere() ;
Vec randomSampleOnHemisphere( Vec const & upDirection ) ;
Vec cosineSampleHemisphere(Vec const &upDirection);


struct Ray {
    Vec o, d;
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {
        d.normalize();
    }
};
struct RayDifferentials : public Ray{
    bool hasDifferentials;
    Vec ox, oy, dir_x, dir_y;
    RayDifferentials(Vec o_, Vec d_, Vec ox_, Vec oy_, Vec dir_x_, Vec dir_y_):
    Ray(o_, d_),
    ox(ox_),
    oy(oy_),
    dir_x(dir_x_),
    dir_y(dir_y_){
        d.normalize();
        dir_x.normalize();
        dir_y.normalize();
        hasDifferentials = true;
    }
    RayDifferentials(Vec o_, Vec d_) : Ray(o_, d_) {
        d.normalize();
        hasDifferentials = false;
    }
};

struct Vec2 {
    double x, y;
    Vec2(double x_=0, double y_=0){ x=x_; y=y_;}
    Vec2 operator+(const Vec2 &b) const { return Vec2(x+b.x,y+b.y); }
    Vec2 operator-(const Vec2 &b) const { return Vec2(x-b.x,y-b.y); }
    Vec2 operator*(double b) const { return Vec2(x*b,y*b); }
    double dot(const Vec2 &b) const { return x*b.x+y*b.y; }
};
template<typename T> Vec2 operator*(T b, Vec2& v){
    return Vec2(v.x*b,v.y*b);
}

struct Mat2{
    double a,b,c,d;
    Mat2(double a_ =0, double b_=0, double c_=0, double d_=0):a(a_), b(b_), c(c_), d(d_){}
    Mat2 operator+(const Mat2 &M) const{return Mat2(a+M.a, b+M.b, c+M.c, d+M.d);}
    Mat2 operator-(const Mat2 &M) const{return Mat2(a-M.a, b-M.b, c-M.c, d-M.d);}
    Mat2 operator*(double s) const{return Mat2(s*a, s*b, s*c, s*d);}
    double& operator[](int i){
        if(i==0) return a;
        else if(i==1) return b;
        else if(i==2) return c;
        else return d;

    }
    double operator[](int i) const{
        if(i==0) return a;
        else if(i==1) return b;
        else if(i==2) return c;
        else return d;

    }
    Mat2 mult(const Mat2 &M) const{
        return Mat2(a*M.a +b*M.c, a*M.b+ b*M.d,
                    c*M.a + d*M.c, c*M.b + d*M.d );
        
    }
    Mat2 inv(){
        return  Mat2(d, -b, -c , a)*(1./(a*d - b*c));
    }
    
    Mat2 T(){
        return Mat2(a,c,b,d);
    }
    double det() const{
        return a*d - b*c;
    }
    
    Vec2 mult(const Vec2& v){
        return Vec2(a*v.x + b*v.y, c*v.x + d*v.y);
    }
    
};

//Row-major storage
struct Mat3{
    double data[9];
    Mat3(double* data_){
        for(int i = 0; i<9; i++){
            data[i] = data_[i];
        }
    }
    double& operator[](int id){
        return data[id];
    }
    Mat3 operator+(const Mat3 &M) const{
        double datanew[9];
        for(int i = 0; i<9; i++){
            datanew[i] = data[i] + M.data[i];
        }
        return Mat3(datanew);
        
    }
    Mat3 operator-(const Mat3 &M) const{
        double datanew[9];
        for(int i = 0; i<9; i++){
            datanew[i] = data[i] - M.data[i];
        }
        return Mat3(datanew);
        
    }
    Mat3 operator*(double s) const{
        double datanew[9];
        for(int i = 0; i<9; i++){
            datanew[i] = data[i]*s;
        }
        return Mat3(datanew);
        
    }
    Mat3 mult(const Mat3 &M) const{
        double datanew[9] = {0,0,0,0,0,0,0,0,0};
        for(int i = 0; i<3; i++){
            for(int j = 0; j<3; j++){
                for (int k = 0; k<3; k++){
                    datanew[3*i + j] += data[3*i+ k]*M.data[3*k +j];
                }
            }
        }
        return Mat3(datanew);
    }
    
    
    Vec mult(const Vec& v){
        float datavec[3] = {0,0,0};
        for(int i = 0; i<3; i++){
            for (int k = 0; k<3; k++){
                datavec[i] += data[3*i+ k]*v[i];
            }
        }
        return Vec(datavec[0], datavec[1], datavec[2]);
    }
    
    static Mat3 Id(){
        double data[9] = {1,0,0,0,1,0,0,0,1};
        return Mat3(data);
    }
    
//    static Mat3 rotate(Vec axis, double angle){
//        double data[9] ={0,0,0,0,0,0,0,0,0};
//        
//    }
    
};





#endif /* RayVec_hpp */
