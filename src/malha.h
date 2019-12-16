#ifndef MALHAH
#define MALHAH
#include "hitable.h"
#include "vec2.h"
#include <vector>
#include <fstream>
#include <math.h>
#include <cstring>
#include <iostream>
#include <sstream>

struct triangle {
public:
    vec3 vert[3];
    vec2 uv[3];
    vec3 normal[3];

    triangle( const std::vector <vec3> &v, const std::vector<vec2> &t, const std::vector<vec3> &n)
    {
        vert[0] = v[0];
        vert[1] = v[1];
        vert[2] = v[2];
        uv[0] = t[0];
        uv[1] = t[1];
        uv[2] = t[2];
        normal[0] = n[0];
        normal[1] = n[1];
        normal[2] = n[2];
    }

    ~triangle(){}
};

class malha: public hitable  {
    public:
        malha() {}
        malha(material *m, const char* path) : mat_ptr(m)  {
            read_OBJ(path);
            // std::cout << "cheguei aqui"<<std::endl;
        }
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        std::vector< triangle > triangles;
        material *mat_ptr;
        bool read_OBJ(const char* path) 
	{

		std::vector< unsigned int > vertexIndices, uvIndices, normalIndices;
		std::vector< vec3 > temp_vertices;
		std::vector< vec2 > temp_uvs;
		std::vector< vec3 > temp_normals;

		std::ifstream f(path);
		if (!f.is_open())
		{
			// std::cout << "File cannot be oppened or does not exist\n";
			return false;
		}

		// std::cout << "file was  oppened!\n";

		
		while (!f.eof())
		{
			char line[256];
			f.getline(line, 256);

			std::stringstream s;
			s << line;

			char junk;

			if ( line[0] == 'v')
			{
				if (line[1] == 't') 
				{
					vec2 uv;
					s >> junk >> junk >> uv[0] >> uv[1];
					temp_uvs.push_back(uv);
				}
				if (line[1] == 'n') 
				{
					vec3 normal;
					s >> junk >> junk >> normal[0] >> normal[1] >> normal[2];
					temp_normals.push_back(normal);
				}
				else {
					vec3 vertex;
					s >> junk >> vertex[0] >> vertex[1] >> vertex[2];
					temp_vertices.push_back(vertex);
				}
			}

			else if ( line[0] == 'f')
			{
				std::string vertex1, vertex2, vertex3;
				unsigned int vertexIndex[3], uvIndex[3], normalIndex[3];

				s >> junk >> vertex1 >> vertex2 >> vertex3;
				int fstslash = vertex1.find("/");
				int sndslash = vertex1.find("/", fstslash+1);
				int trdslash = vertex1.find("/", sndslash+1);
				std::string fst = vertex1.substr(0, fstslash);
				std::string snd = vertex1.substr(fstslash+1, sndslash - fstslash - 1);
				std::string trd = vertex1.substr(sndslash+1);
				vertexIndex[0] = atoi( fst.c_str() );
				uvIndex[0] = atoi(snd.c_str());
				normalIndex[0] = atoi( trd.c_str() );

				fstslash = vertex2.find("/");
				sndslash = vertex2.find("/", fstslash + 1);
				trdslash = vertex2.find("/", sndslash + 1);
				fst = vertex2.substr(0, fstslash);
				snd = vertex2.substr(fstslash + 1, sndslash - fstslash - 1);
				trd = vertex2.substr(sndslash + 1);
				vertexIndex[1] = atoi(fst.c_str());
				uvIndex[1] = atoi(snd.c_str());
				normalIndex[1] = atoi(trd.c_str());

				fstslash = vertex3.find("/");
				sndslash = vertex3.find("/", fstslash + 1);
				trdslash = vertex3.find("/", sndslash + 1);
				fst = vertex3.substr(0, fstslash);
				snd = vertex3.substr(fstslash + 1, sndslash - fstslash - 1);
				trd = vertex3.substr(sndslash + 1);
				vertexIndex[2] = atoi(fst.c_str());
				uvIndex[2] = atoi(snd.c_str());
				normalIndex[2] = atoi(trd.c_str());

				vertexIndices.push_back(vertexIndex[0]);
				vertexIndices.push_back(vertexIndex[1]);
				vertexIndices.push_back(vertexIndex[2]);
				uvIndices.push_back(uvIndex[0]);
				uvIndices.push_back(uvIndex[1]);
				uvIndices.push_back(uvIndex[2]);
				normalIndices.push_back(normalIndex[0]);
				normalIndices.push_back(normalIndex[1]);
				normalIndices.push_back(normalIndex[2]);
			}
		}
		
		for (unsigned int i = 0; i < vertexIndices.size(); i+=3)
		{
			unsigned int v1 = vertexIndices[i];
			unsigned int v2 = vertexIndices[i+1];
			unsigned int v3 = vertexIndices[i+2];

			unsigned int n1 = normalIndices[i];
			unsigned int n2 = normalIndices[i+1];
			unsigned int n3 = normalIndices[i+2];

			unsigned int uv1 = uvIndices[i];
			unsigned int uv2 = uvIndices[i+1];
			unsigned int uv3 = uvIndices[i+2];

			std::vector<vec3> vertices;
			vertices.push_back(temp_vertices[v1 - 1]);
			vertices.push_back(temp_vertices[v2 - 1]);
			vertices.push_back(temp_vertices[v3 - 1]);

			std::vector<vec2> uvs;
			if( uvIndices.size() > 0 ){
				uvs.push_back(temp_uvs[uv1 - 1]);
				uvs.push_back(temp_uvs[uv2 - 1]);
				uvs.push_back(temp_uvs[uv3 - 1]);
			}

			std::vector<vec3> normals;
			if( normalIndices.size() > 0){

				normals.push_back(temp_normals[n1 - 1]);
				normals.push_back(temp_normals[n2 - 1]);
				normals.push_back(temp_normals[n3 - 1]);
			}

			triangle t(vertices, uvs, normals);
			triangles.push_back(t);
		}

		std::cout << "vertSize = " << vertexIndices.size() << "\n";
		std::cout << "normalSize = " << normalIndices.size() << "\n";
		std::cout << "uvSize = " << uvIndices.size() << "\n";

		return true;
	}
    

    bool intersect_triangle(const ray& r, const triangle &tr, float t_min, float t_max, hit_record& rec)const{
        vec3 vertex0 = tr.vert[0];
        vec3 vertex1 = tr.vert[1];  
        vec3 vertex2 = tr.vert[2];
        vec3 edge1, edge2, h, s, q;
        float a,f,u,v;
        edge1 = vertex1 - vertex0;
        edge2 = vertex2 - vertex0;
        h = cross(r.direction(), edge2);
        a = dot( edge1, h );
        if (a > -0.000001 && a < 0.000001){
            return false;    // This ray is parallel to this triangle.
        }
        f = 1.0/a;
        s = r.origin() - vertex0;
        u = f *dot(s,h);
        if (u < 0.0 || u > 1.0)
            return false;
        q = cross(s, edge1);
        v = f * dot(r.direction(), q);
        if (v < 0.0 || u + v > 1.0)
            return false;
        // At this stage we can compute t to find out where the intersection point is on the line.
        float t = f * dot(edge2, q);
        if (t >= t_min && t <= t_max) // ray intersection
        {
            float w = 1.0 - u - v;
            vec3 normal = w*tr.normal[0] + u*tr.normal[1] + v*tr.normal[2]; 
            rec.p = r.origin() + r.direction() * t;
            rec.t = t;
            rec.normal = normal;
            rec.mat_ptr = this->mat_ptr;
            return true;
        }
        else {
            return false;
        }// This means that there is a line intersection but not a ray intersection.
        // compute plane's normal

    }

};


bool malha::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
    bool intersect = false;
    hit_record temp;
    for (int i = 0; i < triangles.size(); ++i)
    {
        bool a = intersect_triangle(r, triangles[i], t_min, t_max, temp);
        if ( a == true)
        {
            
            if (rec.t > temp.t){
                rec = temp;
            }
            intersect = true;
            
        }
    }
    return intersect;
}


#endif
