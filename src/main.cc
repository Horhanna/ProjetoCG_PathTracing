#include <iostream>
#include "sphere.h"
#include "hitable_list.h"
#include "float.h"
#include "camera.h"
#include "material.h"
#include "random.h"
#include "malha.h"
#include <fstream>
#include <algorithm>

vec3 color(const ray& r, hitable *world, int depth) {
    hit_record rec;
    if (world->hit(r, 0.001, 1000, rec)) {
        ray scattered;
        vec3 attenuation;
        vec3 emitted = rec.mat_ptr->emitted();
        if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
             return emitted+attenuation*color(scattered, world, depth+1);
        }
        else {
                return emitted;
            return vec3 (0,0,0);
        }
    }
    //background
    else {
            return vec3(0,0,0);
        vec3 unit_direction = unit_vector(r.direction());
        float t = 0.5*(unit_direction.y() + 1.0);
        return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
    }
}


hitable *random_scene() {
       int n = 8;
       
    hitable **list = new hitable*[n];

    list [0] = new malha( new lambertian(vec3(0.73, 0.73, 0.73)), "./objs/back.obj", "./textura.png");

    list [1] = new malha( new lambertian(vec3(0.73, 0.73, 0.73)), "./objs/box-only.obj");

    list [2] = new malha( new lambertian(vec3(0.73, 0.73, 0.73)), "./objs/ceiling.obj");

    list [3] = new malha( new lambertian(vec3(0.73, 0.73, 0.73)), "./objs/ground.obj");

    list [4] = new malha( new lambertian(vec3(0.65, 0.05, 0.05)), "./objs/left_wall.obj");

    list [5] = new malha( new light(vec3(1, 1, 1), 5), "./objs/light.obj");

    list [6] = new malha( new lambertian(vec3(0.05, 0.65, 0.05)), "./objs/right_wall.obj");

    list [7] = new sphere(vec3(0,0.3, 0), 0.3, new lambertian(vec3(0.05, 0.65, 0.05)));

    return new hitable_list(list,n);
}


int main() {
    int nx = 320;
    int ny = 320;
    int ns = 600;
    std::cout << "P3\n" << nx << " " << ny << "\n255\n";
    hitable *world = random_scene();

    vec3 lookfrom(0,1,3.2);
    vec3 lookat(0,1,1);
    float dist_to_focus = 10.0;
    float aperture = 0.1;
    
//origem= lookfrom, lookat = ponto q a camera ta olhando,fov = 20 ,AR
    camera cam(lookfrom, lookat, vec3(0,-1,0), 50, float(nx)/float(ny), aperture, dist_to_focus);
    std:: vector <vec3> buffer;
    buffer.reserve(nx*ny);
    for (int j = ny-1; j >= 0; j--) {
        for (int i = 0; i < nx; i++) {
            vec3 col(0, 0, 0);
            for (int s=0; s < ns; s++) {
                float u = float(i + random_double()) / float(nx);
                float v = float(j + random_double()) / float(ny);
                ray r = cam.get_ray(u, v);
                col += color(r, world,0);
            }
            col /= float(ns);
            col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
            buffer[j*nx + i] = col;
        }
    }
    std::ofstream ofs3("./imagem_3.ppm", std::ios::out | std::ios::binary); 
    ofs3 << "P6\n" << nx << " " << ny << "\n255\n"; 
    for (unsigned i = 0; i < nx * ny; ++i)
    { 
        
        ofs3 << (unsigned char)(std::min(float(1), (float)buffer[i].x()) * 255) << 
                (unsigned char)(std::min(float(1), (float)buffer[i].y()) * 255) << 
                (unsigned char)(std::min(float(1), (float)buffer[i].z()) * 255); 
    }
    ofs3.close();
}
