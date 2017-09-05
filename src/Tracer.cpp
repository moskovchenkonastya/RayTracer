#include "Tracer.h"
#include "Scene.h"
#include <cstdio>
#include <iostream>
#include <cmath>
#include "glm/gtx/norm.hpp"
#include "glm/gtx/perpendicular.hpp"

#define G 6.674e-11 
#define c 3e+8
#define pi 3.14

using namespace glm;
using namespace std;

SRay CTracer::MakeRay(glm::uvec2 pixelPos, int l1, int l2)
{
  // генерация луча
  SRay ray;

  int W = m_camera.m_resolution.x;
  int H = m_camera.m_resolution.y;

  m_camera.m_viewAngle = vec2(radians(60.0), radians(60.0));
  
  vec3 v_right	= normalize(m_camera.m_right) * (float)((pixelPos.x + l1) / W - 0.5) * l2Norm(m_camera.m_forward) * tan(m_camera.m_viewAngle.x / 2);
  vec3 v_up		= normalize(m_camera.m_up) * (float)((pixelPos.y + l2) / H - 0.5) * l2Norm(m_camera.m_forward) * tan(m_camera.m_viewAngle.y / 2);
   
  ray.m_start = m_camera.m_pos;
  ray.m_dir =  normalize(m_camera.m_forward + v_right + v_up);
  return  ray;
}

float plane_intersetion(SRay ray, double R_disk)
{	
	if (ray.m_start.z * (ray.m_start.z + ray.m_dir.z) > 0)
		return -1.0f;

	double t = -ray.m_start.z / ray.m_dir.z;
	
	float X	= ray.m_start.x +  (float)t * ray.m_dir.x;
	float Y	= ray.m_start.y +  (float)t * ray.m_dir.y;
	return (X * X + Y * Y <= R_disk * R_disk) ? t : -1.0f;
}

float sphere_intersection(SRay ray, double R)
{
	double a = dot(ray.m_dir, ray.m_dir);
	double b = dot(ray.m_start, ray.m_dir);
	double f = dot(ray.m_start, ray.m_start) - R * R;
	double d = (b * b) - (a *  f);
	//cout<< d << endl;

	if(d >= 0)
		return (-b - sqrt(d)) / a <= 1? (-b - sqrt(d)) / a: -1.0f;
	else
		return -1.0f;

}

glm::vec3 CTracer::get_color_plane(SRay ray, double R_disk, float t)
{	
	vec3 color(0,0,0);
	
	vec3 P	= ray.m_start +  (float)t * ray.m_dir;
	if((P.x >= -R_disk && P.x <= R_disk) && (P.y >= -R_disk && P.y <= R_disk))
	{
		int i = (P.y + R_disk) * m_pScene->disk->GetHeight() / (2 * R_disk);
		int j = (P.x + R_disk) * m_pScene->disk->GetWidth() / (2 * R_disk);
		auto pData = (unsigned char*)m_pScene->disk->GetBits();
		int pitch = m_pScene->disk->GetPitch();
		unsigned char b = pData[i * pitch + j * 4];
		unsigned char g = pData[i * pitch + j * 4 + 1];
		unsigned char r = pData[i * pitch + j * 4 + 2];
		unsigned char alpha = pData[i * pitch + j * 4 + 3];
		if(alpha != 0){
			color = vec3(r/255.0, g/255.0, b/255.0);
		}else
			color = vec3(0, 0, 0);
	}
		return color;
	
}

glm::vec3 CTracer::get_color_background(SRay ray)
{	
	vec3 color(0,0,0);
	ray.m_dir = normalize(ray.m_dir);
	int i = (asin(ray.m_dir.z) + pi / 2) * m_pScene->sphere->GetHeight() /  pi;
	int j = (atan2(ray.m_dir.x, ray.m_dir.y) + pi) * m_pScene->sphere->GetWidth() / (2 * pi);

	auto pData = (unsigned char*)m_pScene->sphere->GetBits();
	int pitch = m_pScene->sphere->GetPitch();
	unsigned char b = pData[i * pitch + j * 3];
	unsigned char g = pData[i * pitch + j * 3 + 1];
	unsigned char r = pData[i * pitch + j * 3 + 2];
	
	color = vec3(pow((float)(r / 255.0) , 1.5f), pow((float)(g / 255.0), 1.5f), pow((float)(b / 255.0), 1.5f));
	//color = vec3((float)r / 255.0, (float)g/ 255.0, (float)b / 255.0);
	return color;
	
}



glm::vec3 CTracer::TraceRay(SRay ray, double R, double k, double M)
{	
	vec3 color(0, 0, 0);

	double R_disk = k * R;
	float t = 0;

	vec3 aks;
	float delta_t = 10;
	vec3 v0 = ray.m_dir;
	vec3 pos = ray.m_start;


	float t1, t2;
	for(int i = 0; i < 1000; i++)
	{	
		ray.m_dir = v0;
		ray.m_start = pos;

		aks = -normalize(ray.m_start) * float(G * M / pow(l2Norm(ray.m_start), 2));
		aks = perp(aks, ray.m_dir);
		
		pos = ray.m_start + ray.m_dir * delta_t * float(c) + aks * pow(delta_t, 2) / 2.0f;
	
		v0 = normalize(ray.m_dir * float(c)  + aks * delta_t);

		
		SRay tmp_ray = ray;
		
		ray.m_dir = pos - ray.m_start;

		t1 = plane_intersetion(ray, R_disk);

		t2 = sphere_intersection(ray, R);	

		ray = tmp_ray;
		if ((t1 > 0) && (t2 > 0))
		{

			if (t1 < t2)
			{
				vec3 P	= ray.m_start +  float(t1) * (pos - ray.m_start); 
				if (l2Norm(ray.m_start - P) < l2Norm(ray.m_start - pos)) {
					color = get_color_plane(ray, R_disk, t1);
					return color;
				}
			}
			else
			{
				vec3 P	= ray.m_start + float(t2) * (pos - ray.m_start);
				if (l2Norm(ray.m_start - P) < l2Norm(ray.m_start - pos)) return vec3(0, 0, 0);
			}
		}

		vec3 P;
		if(t2 > 0.0f){
			vec3 P	= ray.m_start + (float)t2 * (pos - ray.m_start);
			if (l2Norm(ray.m_start - P) < l2Norm(ray.m_start - pos)) return vec3(0, 0, 0); 
		}

		if(t1 > 0.0f){
			vec3 P	= ray.m_start +  (float)t1 * (pos - ray.m_start); 
			if (l2Norm(ray.m_start - P) < l2Norm(ray.m_start - pos)) {
				color = get_color_plane(ray, R_disk, t1);
				return color; 
			}		
		}
			
		
	}
	color = get_color_background(ray);
	return color;
  // обратная трассировка  
}

void CTracer::RenderImage(int xRes, int yRes, double M, double k)
{
/*
  std::vector<glm::vec3> disk = GetPixels("data/disk_24.png");
  std::vector<glm::vec3> background = GetPixels("data/stars.jpg");
*/


  m_pScene->disk = LoadImageFromFile("data/disk_32.png");
   m_pScene->sphere = LoadImageFromFile("data/stars.jpg");

  if(m_pScene->disk->GetBPP() == 32)
  {
    auto pData = (unsigned char*)m_pScene->disk->GetBits();
    auto pCurrentLine = pData;
    int pitch = m_pScene->disk->GetPitch();

    for(int i = 0; i < m_pScene->disk->GetHeight(); i++) // Image lines
    {
      for(int j = 0; j < m_pScene->disk->GetWidth(); j++) // Pixels in line
      {	
		
			unsigned char b = pCurrentLine[i * pitch + j * 4];
			unsigned char g = pCurrentLine[i * pitch + j * 4 + 1];
			unsigned char r = pCurrentLine[i * pitch + j * 4 + 2];
			unsigned char alpha = pCurrentLine[i * pitch + j * 4 + 3];
		
        //unsigned char alpha = pCurrentLine[j * pitch + i * 4 + 3];
      }
	  
    }

		
  }
  // расчитаем радиус 
  double R = 2 * G * M / (c * c);
  double R_disk = k * R;
  double z =  5 * R_disk;

  // Rendering
  m_camera.m_resolution = uvec2(xRes, yRes);
  m_camera.m_pixels.resize(xRes * yRes);
  
  // задать параметры камеры прямо
  /*
  m_camera.m_forward	= m_camera.m_up * m_camera.m_right;
  m_camera.m_pos		= vec3(0.0, 1e+11, 2e+9);
  m_camera.m_up			= vec3(0.0, 0.0, 2.0);
  m_camera.m_right		= vec3(3.0, 0.0, 0.0);
   
  //m_camera.m_pos = m_camera.m_forward * (float)(-z);
  m_camera.m_right = vec3(0, -1, 0);
  m_camera.m_forward = vec3(-4, 0, -1);
  m_camera.m_pos = vec3(z, 0, z/4);
  m_camera.m_up = vec3(-1, 0, 4);
  */
  //m_camera.m_up = normalize(m_camera.m_up);
 /*
  m_camera.m_forward	= vec3(-1, 1, -1);
  m_camera.m_pos		= m_camera.m_forward * (float)(-z); // vec3(-12e+10, 0, 2e+10);
  m_camera.m_up			= vec3(1, 1, 0);
  m_camera.m_right		= vec3(1.0, -1, -2);
  
  //m_pScene->norma		= vec3(1.0, -.0, 0.0);
  
  vec3 color(0,0,0);
  

  m_camera.m_forward	= vec3(0.0, 0.0, 1.0);
  m_camera.m_pos		= vec3(0.0, 0.0, -z);
  m_camera.m_up			= vec3(1.0, 0.0, 0.0);
  m_camera.m_right		= vec3(0.0, -1.0, 0.0);

  */
  m_camera.m_forward = vec3(-0.71, 0.69, -0.14); 
  m_camera.m_pos = m_camera.m_forward * (float)(-z); 
  m_camera.m_up = vec3(-0.1, 0.1, 0.99);
  m_camera.m_right = vec3(0.7, 0.72, 0.0);
  SRay ray;
  int l1, l2;
  vec3 color(0, 0, 0);
#pragma omp parallel for
  for(int i = 0; i < yRes; i++) {
	  cout << i << " line complete\n\r";
    for(int j = 0; j < xRes; j++)
    {

	  m_camera.m_pixels[i * xRes + j] = vec3(0);	
	  ray = MakeRay(uvec2(j, i), 0.5, 0.5);
      color =  TraceRay(ray, R, k, M);
	  m_camera.m_pixels[i * xRes + j] += color; 
	  ray = MakeRay(uvec2(j, i), 0.25, 0.25);
	  color =  TraceRay(ray, R, k, M);
	  m_camera.m_pixels[i * xRes + j] += color; 
	  ray = MakeRay(uvec2(j, i), 0.75, 0.75);
      color =  TraceRay(ray, R, k, M);
	  m_camera.m_pixels[i * xRes + j] += color; 
	  ray = MakeRay(uvec2(j, i), 0.25, 0.75);
      color =  TraceRay(ray, R, k, M);
	  m_camera.m_pixels[i * xRes + j] += color; 
	  ray = MakeRay(uvec2(j, i), 0.75, 0.25);
      color =  TraceRay(ray, R, k, M);
	  m_camera.m_pixels[i * xRes + j] += color; 

	  
	 
    }
  }
}


void CTracer::SaveImageToFile(std::string fileName)
{
  CImage image;

  int width = m_camera.m_resolution.x;
  int height = m_camera.m_resolution.y;

  image.Create(width, height, 24);
    
	int pitch = image.GetPitch();
	unsigned char* imageBuffer = (unsigned char*)image.GetBits();

	if (pitch < 0)
	{
		imageBuffer += pitch * (height - 1);
		pitch =- pitch;
	}

	int i, j;
	int imageDisplacement = 0;
	int textureDisplacement = 0;

	for (i = 0; i < height; i++)
	{
    for (j = 0; j < width; j++)
    {
      vec3 color = m_camera.m_pixels[textureDisplacement + j];

      imageBuffer[imageDisplacement + j * 3] = clamp(color.b, 0.0f, 1.0f) * 255.0f;
      imageBuffer[imageDisplacement + j * 3 + 1] = clamp(color.g, 0.0f, 1.0f) * 255.0f;
      imageBuffer[imageDisplacement + j * 3 + 2] = clamp(color.r, 0.0f, 1.0f) * 255.0f;
    }

		imageDisplacement += pitch;
		textureDisplacement += width;
	}

  image.Save(fileName.c_str());
	image.Destroy();
}

CImage* CTracer::LoadImageFromFile(std::string fileName)
{
  CImage* pImage = new CImage;

  if(SUCCEEDED(pImage->Load(fileName.c_str())))
    return pImage;
  else
  {
    delete pImage;
    return NULL;
  }
}