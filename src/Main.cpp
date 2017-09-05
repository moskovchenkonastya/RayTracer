#include "Tracer.h"
#include <stdio.h>


void main(int argc, char** argv)
{
  CTracer tracer;
  CScene scene;

  int xRes = 512;  // Default resolution
  int yRes = 512;
  double M = 8.57e+36;
  double k = 5;

  if(argc == 4) // There is input file in parameters
  {
    FILE* file = fopen(argv[1], "r");
    if(file)
    {
      int xResFromFile = 0;
      int yResFromFile = 0;
	  long double MFromFile = 0;
	  int kFromFile = 0;
      if(fscanf(file, "%d %d %llf %d", &xResFromFile, &yResFromFile, &MFromFile, &kFromFile) == 4)
      {
        xRes = xResFromFile;
        yRes = yResFromFile;
		M = MFromFile;
		k = kFromFile;
      }
      else
        printf("Invalid config format! Using default parameters.\r\n");

      fclose(file);
    }
    else
      printf("Invalid config path! Using default parameters.\r\n");
  }
  else
    printf("No config! Using default parameters.\r\n");

  tracer.m_pScene = &scene;
  tracer.RenderImage(xRes, yRes, M, k);
  tracer.SaveImageToFile("Result.png");
}