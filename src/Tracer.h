#pragma once

#include "glm/glm.hpp"
#include "glm/gtx/norm.hpp"
#include "Types.h"
#include "Scene.h"
#include "omp.h"

#include "string"
#include "atlimage.h"

class CTracer
{
public:
  SRay MakeRay(glm::uvec2 pixelPos, int l1, int l2);  // Create ray for specified pixel
  glm::vec3 TraceRay(SRay ray, double R, double k, double  M); // Trace ray, compute its color
  void RenderImage(int xRes, int yRes, double M, double k);
  void SaveImageToFile(std::string fileName);
  CImage* LoadImageFromFile(std::string fileName);
  glm::vec3 get_color_plane(SRay ray, double R_disk, float t);
  glm::vec3 CTracer::get_color_background(SRay ray);
  
 

public:
  SCamera m_camera;
  CScene* m_pScene;
  Sphere sphere;
};