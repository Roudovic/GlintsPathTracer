//
//  brdf.cpp
//  PathTracer
//
//  Created by Ludovic Theobald on 06/11/2019.
//  Copyright © 2019 Ludovic Theobald. All rights reserved.
//

#include "brdf.hpp"
#include <array>
#include <math.h>


float GGX::G(Vec wi, Vec wo, const Intersection& isect){
    Vec n = isect.get_n();
    float G = G1(wi,n)*G1(wo,n);
//    std::cerr << G << std::endl;
    return G;
}
float GGX::G1(Vec w,  Vec n ){
    float nw = n.dot(w);
    return 2*nw/(nw + sqrt(rough*rough + (1 - rough*rough)*nw*nw));
}

float GGX::D(Vec wi, Vec wo, const Intersection& isect){
    Vec n = isect.get_n();
    Vec wh = (wi+wo).normalize();
    float nwh = n.dot(wh);
    
    return (rough*rough) / (M_PI* pow(( 1 + (rough*rough - 1)*nwh*nwh), 2));
}

float GGX::f(Vec wi, Vec wo, const Intersection& isect){
    Vec n = isect.get_n();
    return (D(wi,wo,isect)*F(wi,wo,isect)*G(wi,wo,isect)) / ( 4 * n.dot(wi) * n.dot(wo));
}

Vec GGX::sample_w(Vec wi,
                  Vec n, Vec t){
    double U1 = (double)(rand()) / RAND_MAX;
    double U2 = (double)(rand()) / RAND_MAX;
    double alpha = rough;
    Vec z = n;
    Vec x = t.normalize();
    Vec y = n.cross(t);
    
    //Vec V = Vec(alpha*wi.dot(x), alpha*wi.dot(y), wi.dot(z)).normalize();
    Vec V = (alpha*wi.dot(x)*x + alpha*wi.dot(y)*y + wi.dot(z)*z).normalize();
    
    Vec T1 = (V.dot(z) < 0.9999) ? V.cross(z).normalize() : x;
    Vec T2 = V.cross(T1);
    float a = 1.0 / (1.0 + V.dot(z));
    float r = sqrt(U1);
    float phi = (U2<a) ? U2/a * M_PI : M_PI + (U2-a)/(1.0-a) * M_PI;
    float P1 = r*cos(phi);
    float P2 = r*sin(phi)*((U2<a) ? 1.0 : V.dot(z));
    Vec N = P1*T1 + P2*T2 + sqrt(fmax(0.0, 1.0 - P1*P1 - P2*P2))*V;
    N = (alpha*N.dot(x)*x + alpha*N.dot(y)*y + fmax(0.0, N.dot(z))*z).normalize();
//    N = Vec(alpha*N.x, alpha*N.y, fmax(0.0, N.z)).normalize();
    Vec dir = wi*(-1.0) + N*2*(wi.dot(N));
    
    return dir;
}
Vec GGX::sample_w(Vec wi, const Intersection& isect){
    Vec n = isect.get_n();
    Vec t = isect.get_dpdu().normalize();
    
    return sample_w(wi, n, t);
}

Vec PatchGGX::sample_w(Vec wi, const Intersection& isect){
    int depth = isect.get_depth();
    if(depth > 1){
        return GGX::sample_w(wi,isect);
    }
    else{
        Mat2 sigma_p_inv = computeSigma_p_inv(isect);
        if(isnan(sigma_p_inv.det())){
            return GGX::sample_w(wi,isect);
        }
        Vec2 uv = Vec2(isect.get_u(),isect.get_v());
        //Express sigma_p_inv as sigmax sigmay and rho
        Mat2 sigma_p = sigma_p_inv.inv();
        double sigma_x = sqrt(sigma_p.a), sigma_y = sqrt(sigma_p.d), rho = sigma_p.b/(sigma_x*sigma_y);
        
        //Sample the pixel footprint
        Vec2 uv_sampled;
        gsl_ran_bivariate_gaussian(r, sigma_x, sigma_y, rho, &(uv_sampled.x), &(uv_sampled.y));
        uv_sampled = uv_sampled + uv;
        Vec n = microNormal->queryTexture<Vec>(uv_sampled.x, uv_sampled.y).normalize();
        
        //Add noise to the sampled normal according to the intrisic roughness of the microfacets
        double dx = gsl_ran_gaussian(r, intrinsic_rough);
        double dy = gsl_ran_gaussian(r, intrinsic_rough);
        n = n + Vec(dx,0,0);
        n = n + Vec(0,dy,0);
        n.z = sqrt( 1 - n.x*n.x + n.y*n.y);
        Vec nref = isect.get_n(), t = isect.get_dpdu(), b = isect.get_dpdv();
        n = n.TangentSpaceToWorld(nref, t, b);
        Vec dir = (wi*(-1.0) + n*2*(wi.dot(n))).normalize();
        return dir;
    }
}

float PatchGGX::D_patch(Vec wi, Vec wo, const Intersection& isect){
    Mat2 sigma_p_inv = computeSigma_p_inv(isect);
    if(isnan(sigma_p_inv.det()) || sigma_p_inv.det() < 1e-6){
//        std::cerr << "D Patch error" << sigma_p_inv.det()<< std::endl;

        return GGX::D(wi,wo,isect);
    }
    double u = isect.get_u();
    double v = isect.get_v();
    Vec2 mu_p(u,v);
    float su = microNormal->get_scale_u(), sv = microNormal->get_scale_v();
    int width = microNormal->get_width();
    int height = microNormal->get_height();
    float h = 1;   // Could be a parameter
    int Nh_u = int(width*su*(1.0/h));
    int Nh_v = int(height*sv*(1.0/h));
//    int Nh_u = int(width*(1.0/h));
//    int Nh_v = int(height*(1.0/h));
    double sigma_r2_inv = 1.0/(intrinsic_rough*intrinsic_rough);
    double sigma_h2_inv = (8*log(2.0)*width*width*su*su)/(h*h);
    Vec wh = (wi+wo).normalize(), n = isect.get_n(), t = isect.get_dpdu(), b = isect.get_dpdv();
    //Project wh on the tangent space and create the 2D representation of it
    Vec s_3d = wh;
//    s_3d = wh;
    Vec2 s = Vec2(s_3d.x, s_3d.y);
//    std::cerr<< s.x << " " << s.y << std::endl;
    
    double Dp = 0;
    
    std::array<double, 2> BB = computeBB_p(sigma_p_inv, 8);
    int el_u_start = (u - BB[0]) *Nh_u ;
    int el_v_start = (v - BB[1]) *Nh_v ;
    int el_u_end = (u + BB[0]) *Nh_u ;
    int el_v_end = (v + BB[1]) *Nh_v ;
    
//    //Square fixed box representing the pixel footprint
//    int margin = 1;
//    int el_u_start = u *Nh_u - margin;
//    int el_v_start = v *Nh_v - margin;
//    int el_u_end = u *Nh_u + margin + 1;
//    int el_v_end = v *Nh_v + margin + 1;
    for ( int el_v = el_v_start; el_v < el_v_end; el_v++){
        for ( int el_u = el_u_start; el_u <el_u_end; el_u++){
            double ui = double(el_u)/Nh_u;
            double vi = double(el_v)/Nh_v;
            Vec si_3d = microNormal->queryTexture<Vec>(ui, vi);
            si_3d = si_3d.TangentSpaceToWorld(n, t, b).normalize();
            Vec2 si = Vec2(si_3d.x, si_3d.y);
            // Implementation for FLAT ELEMENTS
            Vec2 mu_i(ui,vi);
            Mat2 sigma_i_inv  = Mat2(sigma_h2_inv, 0, 0, sigma_h2_inv);
            Mat2 sigma_inv = sigma_p_inv + sigma_i_inv;
            Mat2 sigma = sigma_inv.inv();
            Vec2 mu = sigma.mult(sigma_p_inv.mult(mu_p) + sigma_i_inv.mult(mu_i));


            double cp = sqrt(sigma_p_inv.det())/(2*M_PI);
            // c = 1/det(Sigma)*4pi^2*h^2
            double ci = (sigma_r2_inv*sigma_h2_inv)/(4.0*M_PI*M_PI); //*Nh_u*Nh_v
            double ci_ = ci * exp(-(0.5)*(s-si).dot(s-si)*sigma_r2_inv);
            double G_mu_i = ci_ * exp(-(0.5)*(mu - mu_i).dot(sigma_i_inv.mult(mu - mu_i)));
            double G_mu_p = cp * exp(-(0.5)*(mu - mu_p).dot(sigma_p_inv.mult(mu - mu_p)));
            double c = G_mu_p * G_mu_i;
            double Dp_add = c * 2*M_PI*sqrt(sigma.det());
            if (isnan(Dp_add)){
                continue;
            }
            Dp += Dp_add;

            
        }
    }
//    std::cerr << "D Patch end " << std::endl;

    return Dp;
}
float PatchGGX::D(Vec wi, Vec wo, const Intersection& isect){
    return GGX::D(wi,wo,isect);
}

float PatchGGX::f(Vec wi, Vec wo, const Intersection& isect){
    int depth = isect.get_depth();
    Vec n = isect.get_n();
    if(depth > 1)
        return (D(wi,wo,isect)*F(wi,wo,isect)*G(wi,wo,isect)) / ( 4 * n.dot(wi) * n.dot(wo));
    else
        return (D_patch(wi,wo,isect)*F(wi,wo,isect)*G(wi,wo,isect)) / ( 4 * n.dot(wi) * n.dot(wo));
}

Mat2 PatchGGX::computeSigma_p_inv(const Intersection& isect){
    //[Herbert 89]

    Mat2 J = isect.get_fp();
    //[Herbert 89]
//    double A = J.c*J.c + J.d*J.d + 1;
//    double C = J.a*J.a + J.b*J.b + 1;
//    double B = -2*(J.a*J.c + J.b*J.d) ;
////    std::cout << A << std::endl;
//    
//    double inv_F = 1/(A*C - B*B*0.25);
//    A *= inv_F;
//    B *= inv_F;
//    C *= inv_F;
//    std::cout << A << std::endl;

    Mat2 J_inv = J.inv();
//    std::cout << J_inv.d << std::endl;
    Mat2 sigma_p_inv = J_inv.T().mult(J_inv);
//    std::cout << sigma_p_inv.d << std::endl;

    return sigma_p_inv;
}

std::array<double,2> PatchGGX::computeBB_p(const Mat2& sigma_p_inv, double bound){
    
    double bound_x = sqrt(bound/sigma_p_inv.a);
    double bound_y = sqrt(bound/sigma_p_inv.d);
//    std::cout << bound_y << std::endl;
    std::array<double,2> BB = {fmin(bound_x,0.005),fmin(bound_y,0.005)};
    if(isnan(BB[0]) || isnan(BB[1])){
//        std::cout << bound_x << std::endl;
        return {1e-3, 1e-3};
    }
//    if(BB[0]>1e-5 ||BB[1]>1e-5)
//        std::cout << BB[0] << " " << BB[1] << std::endl;


    return BB;
}



