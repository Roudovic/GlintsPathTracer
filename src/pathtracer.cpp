#ifndef pathtracer_hpp
#define pathtracer_hpp

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include "Texture.hpp"
#include "RayVec.hpp"
#include "brdf.hpp"
#include "Intersection.hpp"
#include "SrgbTransform.hpp"
#include "sphere.hpp"
#include "triangle.hpp"
#include "mesh.hpp"
#include <chrono>

#define SETUP 1



Vec refract(Vec const &I, Vec &N, double ior){
    Vec Icpy = I;
    Vec Nrefr = N;
    double cosi = N.dot(Icpy);
    double ni = 1, nr = ior;
    if(cosi < 0){
        cosi = -cosi;
    }
    else{
        Nrefr = -1.0*N;
        std::swap(ni,nr);
    }
    double eta = ni/nr;
    double k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? 0 : eta * Icpy + (eta * cosi - sqrtf(k)) * Nrefr;
    
}

double powerHeuristic(int ns, double pdf_s, int nd, double pdf_d){
    double f = ns * pdf_s, g = nd * pdf_d;
    return (f * f) / (f * f + g * g);
}
double balanceHeuristic(int ns, double pdf_s, int nd, double pdf_d){
    double f = ns * pdf_s, g = nd * pdf_d;
    return f  / (f+g);
}
double rootHeuristic(int ns, double pdf_s, int nd, double pdf_d){
    double f = ns * pdf_s, g = nd * pdf_d;
    return sqrt(f)  / (sqrt(f) + sqrt(g));
}



// Scene :
static std::vector< Shape* > shapes;
static std::vector< unsigned int > lights;
static std::string filename_albedo = "Textures/Brushed_gold_Base_Color.png";
static std::string filename_normal = "Textures/Brushed_gold_Normal.png";
static std::shared_ptr<Texture> texAlbedo = std::make_shared<Texture>(filename_albedo, 1,1,1);
static std::shared_ptr<Texture> texNormal = std::make_shared<Texture>(filename_normal, 1,5,1);


// lights is the whole set of emissive objects
static std::shared_ptr<GGX> ggx_glossy1 = std::make_shared<GGX>(0.2, 0.01);
static std::shared_ptr<GGX> ggx_metal1 = std::make_shared<GGX>(0.95, 0.05);


static std::shared_ptr<Texture> microNormal = std::make_shared<Texture>(filename_normal, 1,3,3);
static std::shared_ptr<PatchGGX> p_ggx1 = std::make_shared<PatchGGX>(  1 , 0.05, 0.005, microNormal);
static std::shared_ptr<Phong> phong1 = std::make_shared<Phong>(1, 3);


inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }
inline bool intersectScene(const RayDifferentials &r, double &t, int &id, Intersection &isect){
    double d, inf=t=1e20;
    Intersection isect_tmp;
    for(int i= 0 ;i < shapes.size(); ++i){
        if((d=shapes[i]->intersect(r, isect_tmp, shapes[i]->needRayDiff()))&&d<t){
            t=d;
            id=i;
            isect = isect_tmp;
        }
    }
    return t<inf;
}



Vec radiance(const RayDifferentials &r, int depth, int sampleIt){
    double t;                               // distance to intersection
    int id=0;                               // id of intersected object
    Intersection isect;
    if (!intersectScene(r, t, id, isect)) return Vec(); // if miss, return black
    double u = isect.get_u();
    double v = isect.get_v();
    const Shape &obj = *(shapes[id]);        // the hit object
    
    Vec x=r.o+r.d*t;                  // the hit position
    
    Vec ng = isect.get_ng();  // the GEOMETRIC normal of the object at the hit
    Vec n;
    if (obj.isNormalTextured()){
        n = obj.getNormalTex()->queryTexture<Vec>(u, v);
        n = n.TangentSpaceToWorld(ng, isect.get_dpdu(), isect.get_dpdv());
        isect.set_n(n);
    }
    else{
        n = ng;
    }
    //f=obj.c;                  // the color of the object at the hit
    Vec f = obj.isAlbedoTextured() ? obj.getAlbedoTex()->queryTexture<Vec>(u, v) : obj.getColor();
    int max_depth = 4;
    if (++depth>max_depth) return Vec(); // we limit the number of rebounds in the scene
    isect.set_depth(depth);
    if (obj.getReflectionType() == EMMISSIVE){ // we hit a light
        return obj.getEmission();
    }
    if (obj.getReflectionType() == DIFFUSE){                  // Object with a scattering surface
        Vec rad = 0;
        int n_bounce = 1;
        
        for( unsigned int n_ray = 0; n_ray<n_bounce ; ++n_ray){
            
            // Multiple importance sampling
            // Sample from the brdf
            Vec wi = r.d * (-1.0);
            wi = wi.normalize();
            Vec dir = obj.getBrdf()->sample_w(wi, isect);
            dir = dir.normalize();
            double cosi = clamp(dir.dot(n));

            double pdf = obj.getBrdf()->pdf(wi, dir, isect);
            RayDifferentials ray2(x+0.01*dir, dir);
            double fr = obj.getBrdf()->f(wi, dir, isect);
            Vec rad_plus_brdf = (radiance(ray2, depth, sampleIt) )*cosi*fr*(1./pdf);
            if(dir.dot(ng) < 0){
                rad_plus_brdf = 0;
            }
            
            //Sample one light at a time
            size_t n_light = lights.size();
            int id_sampled_light = rand()%n_light;
            const Shape &light = *(shapes[ lights[id_sampled_light] ]);
            double pdf_li;
            Vec dir_light = light.randomSample(x, pdf_li).normalize();
            int idInter;
            double dist_inter;
            Intersection  isect_bis;
            RayDifferentials ray_next(x+0.01*dir_light, dir_light);
            intersectScene(ray_next, dist_inter, idInter, isect_bis);
            double cosi_light = clamp(dir_light.dot(n));
            double fr_light = 0;
            Vec rad_light = 0;
            if(idInter == lights[id_sampled_light] ){
                fr_light = obj.getBrdf()->f(wi, dir_light, isect);
                rad_light =  fr_light * light.getEmission() \
                * cosi_light * (1./pdf_li);
                
                double pdf_li_brdf;
                intersectScene(ray2, dist_inter, idInter, isect_bis);
                if(idInter ==lights[id_sampled_light] ){
                    pdf_li_brdf = pdf_li;
                }
                else{
                    pdf_li_brdf = 0;
                }
                double pdf_brdf_li = obj.getBrdf()->pdf(wi, dir_light, isect);
                double w_b = powerHeuristic(1, pdf, 1, pdf_li_brdf);
                double w_li = powerHeuristic(1, pdf_li, 1, pdf_brdf_li);
                if(dir_light.dot(ng) < 0){
                    rad_light = 0;
                }
//                double w_b = 0;
//                double w_li = 1;
                
                rad =  w_b * rad_plus_brdf + w_li *rad_light;
                
            }
            else{
                rad =  rad_plus_brdf;
            }
            rad = rad + obj.getEmission();
            
        }
        if(isnan(rad.dot(rad))){
            rad = Vec(0,0,0);
//                            std::cerr << "Rad is NaN" <<std::endl;
        }
        
        return f.mult(rad);
        
        
    }
    else if (obj.getReflectionType() == MIRROR) {           // Ideal SPECULAR reflection
        Vec dir = r.d - 2*(r.d.dot(n))*n;
        dir = dir.normalize();
        RayDifferentials ray2(x+0.01*dir, dir);
        return obj.getEmission() + radiance(ray2,--depth, sampleIt);
    }
    else if (obj.getReflectionType() == GLASS) {           // Ideal  refraction
        double n1 = 1, n2 = 1.7;
        Vec rad;
        Vec dout;
        double cosi1 = r.d.dot(n);
        if(cosi1*cosi1 > 1 - (n2*n2)/(n1*n1)){
            dout = refract(r.d, n, n2);
        }
        else{
            dout=r.d - 2*(r.d.dot(n))*n;
        }
        RayDifferentials ray2(x+0.01*dout,dout);
        rad = rad + radiance(ray2, --depth, sampleIt);
        return rad;
    }
    
    
    return Vec();
}





int main(int argc, char *argv[]){
    bool is_video = false;
    int n_frame = 30;
    
    std::chrono::time_point<std::chrono::steady_clock> start_point = std::chrono::high_resolution_clock::now();
    
#if SETUP==1
    shapes.reserve(9);
    shapes.push_back( new Sphere(1e5, Vec( 1e5+99,40.8,81.6),    Vec(0,0,0)    ,Vec(.75,.25,.25), DIFFUSE)   );//Right
    shapes.push_back( new  Sphere(1e5, Vec(-1e5+1,40.8,81.6),   Vec(0,0,0)    ,Vec(.25,.25,.75), DIFFUSE,false, false, texAlbedo, texNormal)   );//Left
    shapes.push_back( new Sphere(1e5, Vec(50,40.8, 1e5+170),        Vec(1,1,1)*0    ,Vec(.75,.75,.75), DIFFUSE)   );//Front
    shapes.push_back( new Sphere(1e5, Vec(50,40.8,-1e5),    Vec(1,1,1)*0.0  ,Vec(.95,.95,.95)     , DIFFUSE)   );//Back
    shapes.push_back( new Sphere(1e5, Vec(50, 1e5 + 81.6, 81.6),       Vec(0,0,0)    ,Vec(.75,.75,.75), DIFFUSE)   );//Top
    shapes.push_back( new Sphere(1e5, Vec(50,-1e5,81.6),   Vec(0,0,0)    ,Vec(.75,.75,.75), DIFFUSE ,false, false , 0, 0)   );//Bottom
    //    shapes.push_back( new Sphere(3.5,Vec(29,32.5,67),          Vec(0,0,0)    ,Vec(0.75,0.75,0.75)*.999 , DIFFUSE, true, false , texAlbedo, texNormal, ggx1)    );//Mirr
    shapes.push_back( new Sphere(16.5,Vec(73,16.5,88),          Vec(0,0,0)    ,Vec(1,1,1)*.999 , GLASS)     );//Change to Glass
    shapes.push_back( new Sphere(3,  Vec(50,71.6,65)  ,          Vec(10,10,10)*12.8 ,Vec(0,0,0)      , EMMISSIVE) );//Light
//    shapes.push_back( new Sphere(1,  Vec(15,80.6 - 6,90)  ,          Vec(10,10,10)*2.9   ,Vec(0,0,0)      , EMMISSIVE) );//Light
    //    shapes.push_back(new Triangle(Vec(3,36.5,88),Vec(73,36.5,88), Vec(43,66.5,28), Vec(0,0,0), Vec(.95,.95,.75), DIFFUSE,  true, false , texAlbedo, texNormal, p_ggx1));
    //    shapes.push_back(new Mesh("Models/monkey.off",Vec(30,23.5,58), Mat3::Id(),20,Vec(0,0,0)    ,Vec(0.75,0.75,0.75)*.999 , GLASS, true, false , texAlbedo, texNormal, ggx1  ));
    //    shapes.push_back(new Mesh("Models/monkey.off",Vec(73,46.5,88), Mat3::Id(),20,Vec(0,0,0)    ,Vec(0.75,0.75,0.75)*.999 , GLASS  ));
    //    lights.push_back(3);
    shapes.push_back(new Mesh("Models/cat_simple2.off",Vec(30,0.5,60), Mat3::Id(),1.6,Vec(0,0,0)    ,Vec(1,1,1)*.56 , GLASS ,false, false , texAlbedo, texNormal, ggx_glossy1));
//        shapes.push_back( new Sphere(23.5,Vec(29,32.5,67),          Vec(0,0,0)    ,Vec(0.75,0.75,0.75)*.999 , DIFFUSE, false, false , texAlbedo, texNormal, ggx_glossy1)    );//Mirr
//    shapes.push_back(new Mesh("Models/sphere.off",Vec(30,34.5,60), Mat3::Id(),23.5,Vec(0,0,0)    ,Vec(0.75,0.75,0.75)*.999 , GLASS ,false, false , texAlbedo, texNormal, ggx_metal1));

    
    lights.push_back( 7);
//    lights.push_back( 8 );
    
#elif SETUP==2
    spheres.reserve(6);
    
    spheres.push_back( Sphere(1e5, Vec(50,40.8,-1e5-500),    Vec(1,1,1)*0.1  ,Vec(.50,.80,.999)     , DIFFUSE)   );//Back
    //    spheres.push_back( Sphere(1e5, Vec(50, 1e5 + 81.6, 81.6),       Vec(0,0,0)    ,Vec,(.999,.99,.99), DIFFUSE)   );//Top
    spheres.push_back( Sphere(1e5, Vec(50,-1e5,81.6),   Vec(0,0,0)    ,Vec(.75,.75,.75), DIFFUSE, false, false, 0,0,ggx1)   );//Bottom
    spheres.push_back( Sphere(40.5,Vec(50,45.5,67),          Vec(0,0,0)    ,Vec(0.85,0.85,0.85)*.999 , DIFFUSE, true, false, texAlbedo, texNormal, ggx1)    );//Mirr
    spheres.push_back( Sphere(16.5,Vec(50,25.5,108),          Vec(0,0,0)    ,Vec(1,1,1)*.999 , GLASS)     );//Change to Glass
    spheres.push_back( Sphere(10,  Vec(50 ,150.6 - 6,175)  ,          Vec(10,10,10)*6.699   ,Vec(0,0,0)      , EMMISSIVE) );//Light
    spheres.push_back( Sphere(2,  Vec(5,80.6 - 6,90)  ,          Vec(10,10,10)*10.199   ,Vec(0,0,0)      , EMMISSIVE) );//Light
    lights.push_back( 4 );
    lights.push_back( 5 );
#endif
    
    for(unsigned int frame = 0; frame < n_frame; frame++){
        int spp = 40;
        int w=1024, h=768, samps = spp; // # samples
        Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).normalize()); // camera center and direction
        Vec cx=Vec(w*.5135/h), cy=(cx.cross(cam.d)).normalize()*.5135, *pixelsColor=new Vec[w*h];
        
        
        // ray trace:
#pragma omp parallel for schedule(dynamic)
        for (int y=0; y<h; y++) {                       // Loop over image rows
            fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
            for (unsigned short x=0; x<w; x++) {  // Loop cols
                Vec r(0,0,0);
                
                for( unsigned int sampleIt = 0 ; sampleIt < samps ; ++sampleIt ) {
                    double dx = ((double)(rand())/RAND_MAX);
                    double dy = ((double)(rand())/RAND_MAX);
                    Vec d = (cx*( ( x + dx )/w - .5) +
                             cy*( ( y + dy )/h - .5) + cam.d).normalize();
                    //camera ray differentials
                    Vec dir_x = (cx*( ( x + 1 + dx )/w - .5) +
                                 cy*( ( y + dy )/h - .5) + cam.d).normalize();
                    Vec dir_y = (cx*( ( x  + dx )/w - .5) +
                                 cy*( ( y + 1 + dy )/h - .5) + cam.d).normalize();
                    r = r + radiance(RayDifferentials(cam.o+d*140,d,
                                                      cam.o+dir_x*140,dir_x,
                                                      cam.o+dir_y*140,dir_y),0, sampleIt)*(1./samps);
                    
                }
                //#pragma omp critical
                pixelsColor[x + (h-1-y) * w] = pixelsColor[x + (h-1-y) * w] + Vec(clamp(r.x),clamp(r.y),clamp(r.z));
            } // Camera rays are pushed ^^^^^ forward to start in interior
        }
        
        auto endtime_point = std::chrono::high_resolution_clock::now();
        long long start = std::chrono::time_point_cast<std::chrono::seconds>(start_point).time_since_epoch().count();
        long long end = std::chrono::time_point_cast<std::chrono::seconds>(endtime_point).time_since_epoch().count();
        
        std::cout <<  "\n Rendering took " << (end - start)  <<" seconds. \n" ;
        
        // save image:
        
        std::string image_name;
        if(is_video){
            image_name = "frame_" + std::to_string(frame)+".ppm";
        }
        else{
            image_name = "image.ppm";
        }
        FILE *f = fopen(image_name.c_str(), "w");         // Write image to PPM file.
        fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
        for (int i=0; i<w*h; i++)
            fprintf(f,"%d %d %d ", (int)(SrgbTransform::linearToSrgb(pixelsColor[i].x) * 255), (int)(SrgbTransform::linearToSrgb(pixelsColor[i].y) * 255), (int)(SrgbTransform::linearToSrgb(pixelsColor[i].z) * 255));
        fclose(f);
        delete[] pixelsColor;
        pixelsColor = nullptr;
        if(!is_video){
            break;
        }
    }
    return 0;
}


#endif
