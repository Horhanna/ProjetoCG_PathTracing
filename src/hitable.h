#ifndef HITABLEH
#define HITABLEH
#include "ray.h"

class material;

struct hit_record
{
    float t = 10000;
    vec3 p = vec3 (0,0,0);
    vec3 normal = vec3 (0,0,0);
    material * mat_ptr = nullptr;
};

class hitable  {
    public:
        virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const = 0;
};


#endif
