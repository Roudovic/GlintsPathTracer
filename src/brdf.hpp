//
//  brdf.hpp
//  PathTracer
//
//  Created by Ludovic Theobald on 06/11/2019.
//  Copyright © 2019 Ludovic Theobald. All rights reserved.
//

#ifndef brdf_hpp
#define brdf_hpp

#include <stdio.h>
#include "RayVec.hpp"
#include "Texture.hpp"
#include "Intersection.hpp"
#include <math.h>
#include <gsl/gsl_randist.h>


class BRDF{
public:
    virtual float f(Vec wi, Vec wo, const Intersection& isect) = 0;
    virtual Vec sample_w(Vec wi, const Intersection& isect) = 0;
    virtual float pdf(Vec wi,Vec wo, const Intersection& isect) = 0;
    virtual std::string getName() const{
        return name;
    }
private :
    std::string name;
};

class microfacetBRDF : public BRDF{
public:
    microfacetBRDF(float F0_):F0(F0_){}
    virtual float G(Vec wi, Vec wo, const Intersection& isect) = 0;
    virtual float D(Vec wi, Vec wo, const Intersection& isect) = 0;
    virtual float F(Vec wi, Vec wo, const Intersection& isect){
        Vec wh = (wi+wo).normalize();
        //return 1;
        return F0 + (1-F0)*pow(1 - fmax(0,(wi.dot(wh))) , 5);
    }
    
private:
    float F0;
};

class Lambertian : public BRDF{
public:
    Lambertian(float kd_):kd(kd_){}
    float f(Vec wi, Vec wo, const Intersection& isect){
        return kd/M_PI;
    }
    virtual Vec sample_w(Vec wi, const Intersection& isect){
        Vec n = isect.get_n();
        Vec t = isect.get_dpdu().normalize();
        return sample_w(wi,n,t);
    }
    virtual Vec sample_w(Vec wi, Vec n, Vec t ){
//        return randomSampleOnHemisphere(n);
        return cosineSampleHemisphere(n);
    }
    
    virtual float pdf(Vec wi, Vec wo, const Intersection& isect){
        Vec n = isect.get_n();
        return pdf(wi,wo,n);
    }
    virtual float pdf(Vec wi, Vec wo, Vec n){
        double absCosTheta = abs(wo.dot(n));
//        return 1/M_PI;
        return absCosTheta/(M_PI);
    }
private:
    float kd;
};

class Phong : public BRDF{
public:
    Phong(float ks_, float s_): ks(ks_), s(s_){}
    virtual float f(Vec wi, Vec wo, const Intersection& isect){
        Vec n = isect.get_n();
        Vec wh = (wi+wo).normalize();
        return ks*pow((n.dot(wh)),s);
    }
    virtual Vec sample_w(Vec wi, Vec n, Vec t){
        return randomSampleOnHemisphere(n);
    }
    virtual Vec sample_w(Vec wi,const Intersection& isect){
        Vec n = isect.get_n();
        return randomSampleOnHemisphere(n);
    }
    
    virtual float pdf(Vec wi, Vec wo, const Intersection& isect){
//TODO : Faux mais non utilise
        return 1/M_PI;
    }
private:
    float ks, s;
};

class GGX : public microfacetBRDF{
public:
    GGX(float F0_, float rough_):microfacetBRDF(F0_), rough(rough_){}
    virtual float G(Vec wi, Vec wo, const Intersection& isect);
    float G1(Vec w, Vec n);
    virtual float D(Vec wi, Vec wo, const Intersection& isect);
    virtual float f(Vec wi, Vec wo, const Intersection& isect);
    virtual Vec sample_w(Vec wi, const Intersection& isect);
    virtual Vec sample_w(Vec wi, Vec n, Vec t);
    virtual float pdf(Vec wi, Vec wo, const Intersection& isect){
        Vec wh = (wo+wi).normalize();
        Vec n = isect.get_n();
        float D_ = D(wi,wo,isect);
        float G1_ = G1(wi,wh);
        float Dv = (G1_*fmax(0,wi.dot(wh))*D_ )/ (wi.dot(n));
        return Dv/(4*wi.dot(wh));
    
    }
protected:
    float rough, met;
};

// PatchGGX Material is only used for the first bounce of light (except for bounces going through glass). Beyond that, the effect is negligible.
class PatchGGX : public GGX{
public:
    PatchGGX(float F0_, float rough_, float intrinsic_rough_, std::shared_ptr<Texture> microNormal_):GGX(F0_, rough_){
        microNormal = microNormal_;
        intrinsic_rough = intrinsic_rough_;
    }
    virtual float D_patch(Vec wi, Vec wo, const Intersection& isect);
    virtual float D(Vec wi, Vec wo, const Intersection& isect);
    virtual Vec sample_w(Vec wi, const Intersection& isect);
    Mat2 computeSigma_p_inv(const Intersection& isect);
    std::array<double,2> computeBB_p(const Mat2& sigma_p_inv, double bound);
    virtual float pdf(Vec wi, Vec wo, const Intersection& isect){
        int depth = isect.get_depth();
        Vec wh = (wo+wi).normalize();
        float pdf_;
        if(depth>1){
             pdf_ = GGX::pdf(wi,wo,isect);
        }
        else{
            pdf_ = D_patch(wi,wo,isect)/(4*wi.dot(wh));
        }
        return pdf_;
    }
    //TODO : Distinguish between the first and subsequent bounces 
    virtual float f(Vec wi, Vec wo, const Intersection& isect);
    virtual std::string getName() const{
        return "PatchGGX";
    }
private:
    std::shared_ptr<Texture> microNormal;
    gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
    float intrinsic_rough;
    
    
};


#endif /* brdf_hpp */
