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
#include "lodepng.h"

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
        }
		malha(material *m, const char* path, const char *texture_path) : mat_ptr(m)  {
            read_OBJ(path);
            decodeOneStep(texture_path);
        }
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        std::vector< triangle > triangles;
        material *mat_ptr;
		int texture_width, texture_height;
		std::vector<vec3> texture_buffer;
        bool read_OBJ(const char* path) 
	{

		std::vector< unsigned int > vertexIndices, uvIndices, normalIndices;
		std::vector< vec3 > temp_vertices;
		std::vector< vec2 > temp_uvs;
		std::vector< vec3 > temp_normals;

		std::ifstream f(path);
		if (!f.is_open())
		{
			
			return false;
		}
		
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
		vec3 v0= tr.vert[0];
		vec3 v1= tr.vert[1];
		vec3 v2= tr.vert[2];

		vec3 v01= v1 - v0;
		vec3 v02= v2 - v0;

		vec3 normal= cross(v01, v02);
		
		float parallel= dot(r.direction(), normal); 
		if (fabs (parallel) < 0.0001)
			return false;
		//distancia origem a v0
		float d = dot(normal, v0);
		//distancia ponto ao plano
		float t = -(dot(normal, r.origin())-d)/parallel;
		if (t < t_min)
			return false;
		//intersecao do raio com o plano
		vec3 p = r.origin() + r.direction() * t;

		vec3 vp0 = p - v0;
		//vetor perpendicular
		vec3 c0 = cross(v01, vp0);
		//se os vetores estiverem em direcoes opostas retorna falso
		if (dot(normal, c0)< 0)
			return false;

		vec3 v21 = v2 - v1;
		vec3 vp1 = p - v1;
		vec3 c1 = cross(v21, vp1);
		if (dot(normal, c1)< 0)
			return false;

		v02 = v0 - v2;
		vec3 vp2 = p - v2;
		vec3 c2 = cross(v02, vp2);
		if (dot(normal, c2)<0)
			return false;

		rec.t= t;
		rec.p= p;
		rec.mat_ptr = this-> mat_ptr;

		float a = c0.length()/normal.length();
		float b = c2.length()/normal.length();
		float c = c1.length()/normal.length();

		rec.normal = tr.normal[0]*a + tr.normal[1]*b + tr.normal[2]*c;
		if(texture_buffer.size() > 0){
			vec2 st = a*tr.uv[0] + b*tr.uv[1] + c*tr.uv[2];
			int tx = std::floor( st.x()*(float)texture_width);
			int ty = std::floor( (1.0f - st.y())*(float)texture_height);
			int idx = (ty*texture_width + tx);
			if(idx >= 0 && idx < texture_buffer.size())
				rec.mat_ptr = new lambertian(texture_buffer[idx]);
		}
		return true;
		
		
		// vec3 vertex0 = tr.vert[0];
        // vec3 vertex1 = tr.vert[1];  
        // vec3 vertex2 = tr.vert[2];
        // vec3 edge1, edge2, h, s, q;
        // float a,f,u,v;
        // edge1 = vertex1 - vertex0;
        // edge2 = vertex2 - vertex0;
        // h = cross(r.direction(), edge2);
        // a = dot( edge1, h );
        // if (a > -0.000001 && a < 0.000001){
        //     return false;    
        // }
        // f = 1.0/a;
        // s = r.origin() - vertex0;
        // u = f *dot(s,h);
        // if (u < 0.0 || u > 1.0)
        //     return false;
        // q = cross(s, edge1);
        // v = f * dot(r.direction(), q);
        // if (v < 0.0 || u + v > 1.0)
        //     return false;
        
        // float t = f * dot(edge2, q);
        // if (t >= t_min && t <= t_max) 
        // {
        //     float w = 1.0 - u - v;
        //     vec3 normal = w*tr.normal[0] + u*tr.normal[1] + v*tr.normal[2]; 
        //     rec.p = r.origin() + r.direction() * t;
        //     rec.t = t;
        //     rec.normal = normal;
        //     rec.mat_ptr = this->mat_ptr;
		// 	if(texture_buffer.size() > 0){
		// 		vec2 st = w*tr.uv[0] + u*tr.uv[1] + v*tr.uv[2];
		// 		int tx = std::floor( st.x()*(float)texture_width);
		// 		int ty = std::floor( (1.0f - st.y())*(float)texture_height);
		// 		int idx = (ty*texture_width + tx);
		// 		if(idx >= 0 && idx < texture_buffer.size())
		// 			rec.mat_ptr = new lambertian(texture_buffer[idx]);
		// 	}
        //     return true;
        // }
        // else {
        //     return false;
        // }

    }


	void decodeOneStep(const char* filename) {
			std::vector<unsigned char> png;
			std::vector<unsigned char> image; 
			unsigned width, height;
			lodepng::State state; 

			unsigned error = lodepng::load_file(png, filename); 
			if(!error) error = lodepng::decode(image, width, height, state, png);

			
			if(error) std::cout << "decoder error " << error << ": "<< lodepng_error_text(error) << std::endl;

			texture_buffer.reserve( (int)width*(int)height);
			for(int i = 0; i < image.size(); i+=4)
				texture_buffer.push_back( vec3( float(image[i]), float(image[i+1]), float(image[i+2]) )/255.0f );
			
			texture_width = (int)width;
			texture_height = (int)height;
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
