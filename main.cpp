#include <iostream>

#include "Eigen.h"
#include "ImplicitSurface.h"
#include "Volume.h"
#include "MarchingCubes.h"

using namespace std;

int main()
{
	std::string filenameIn = "../../Data/surface_data/normalized.pcb";
	std::string filenameOut = "result.off";

	// implicit surface
	// TODO: you have to switch between these surface types
	Vector3d center(0.5, 0.5, 0.5);
	float radius = 0.4;
	float small_radius = 0.1;
	ImplicitSurface* surface;
	// surface = new Sphere(center, radius);
	// surface = new Torus(center, radius, small_radius);
	// surface = new Hoppe(filenameIn);
	surface = new RBF(filenameIn);

	// fill volume with signed distance values
	unsigned int resolution = 20; 
	double extra_space = 0.1;
	Vector3d min_point(0, 0, 0);
	Vector3d max_point(1, 1, 1);
	min_point = (min_point.array() - extra_space).matrix();
	max_point = (max_point.array() + extra_space).matrix();
	
	Volume vol(min_point, max_point, resolution, resolution, resolution, 1);

	for (unsigned int x = 0; x < vol.getDimX(); x++)
	{
		for (unsigned int y = 0; y < vol.getDimY(); y++)
		{
			for (unsigned int z = 0; z < vol.getDimZ(); z++)
			{
				Vector3d p = vol.pos(x, y, z);
				double val = surface->Eval(p);
				vol.set(x,y,z, val);
			}
		}
	}

	// extract the zero iso-surface using marching cubes
	SimpleMesh mesh;
	for (unsigned int x = 0; x < vol.getDimX() - 1; x++)
	{
		std::cerr << "Marching Cubes on slice " << x << " of " << vol.getDimX() << std::endl;

		for (unsigned int y = 0; y < vol.getDimY() - 1; y++)
		{
			for (unsigned int z = 0; z < vol.getDimZ() - 1; z++)
			{
				ProcessVolumeCell(&vol, x, y, z, 0.00f, &mesh);
			}
		}
	}

	// write mesh to file
	if (!mesh.WriteMesh(filenameOut))
	{
		std::cout << "ERROR: unable to write output file!" << std::endl;
		return -1;
	}

	delete surface;

	return 0;
}
