#pragma once

#include "Types.h"
#include "atlimage.h"

class Sphere
{	
	glm::vec3 center;
    long double radius;
    //Color color
};

class CScene
{
  // Set of meshes
public:
	CImage *disk;
	CImage *sphere;
    glm::vec3 norma;

};