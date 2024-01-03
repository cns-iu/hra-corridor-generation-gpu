// code is from https://github.com/AhmadM-DL/Point-In-Polyhedron
#ifndef PIP_HELPER_H
#define PIP_HELPER_H

#include <stdio.h>
#include <math.h>

//Loading an Obj File Required Libs
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

//Cuda Libs
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
// #include "device_functions.h"
#include "helper_math.h"
//#include "GpuTimer.h"

// Random Number Generator
#include <random>


/**
* Thie flag corresponds to 4 different states of a cell
* empty(11b)(3d) | singualr(10b)(2d) | inside(01b)(1d) | outside(00b)(0d)
*/
#define CELL_STATE_OFFSET 30
#define CELL_STATE ((uint)3<<CELL_STATE_OFFSET)

//#define CELL_STATE_EMPTY ((uint)3<<CELL_STATE_OFFSET)
#define CELL_STATE_SINGULAR ((uint)2<<CELL_STATE_OFFSET)
#define CELL_STATE_INSIDE ((uint)1<<CELL_STATE_OFFSET)
#define CELL_STATE_OUTSIDE ((uint)0<<CELL_STATE_OFFSET)

/**
* The triangle list contains a sublist for each occupied voxel.
* This flag is set in each element that terminates a sublist.
*/
#define TERMINATE_SUBLIST ((uint)1<<31)

#define MAX_SEARCH_UN_SINGULAR 2




// A routine to calculate the dot product of two vertices
// v1 and v2 should be of length 3
inline float __host__ __device__  dot_product(float *v1, float *v2) {
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// A routine to calculate the dot product of two vertices
// v1 and v2 should be of length 3
inline float __host__ __device__  dot_product(double *v1, double *v2) {
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// A routine to calculate the cross product of two vertices
// v1 and v2 should be of length 3
inline void __host__ __device__ cross_product(float *v1, float *v2, float *cross) {
	cross[0] = v1[1] * v2[2] - v1[2] * v2[1];
	cross[1] = -1 * (v1[0] * v2[2] - v1[2] * v2[0]);
	cross[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

inline __host__ void print_if_cuda_error(int line) {
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess) {
		printf("\nError On Line %d: %s\n", line - 1, cudaGetErrorString(err));
	}
}


// A routin to laod an object (.obj) file into an array of triangle, vertices, numv (number of vertices), and numtri (number of triangles)
void objLoader(const char *path, std::vector<int> *triangles_vector, std::vector<float> *vertices_vector, int* numv, int* numtri, float3* min, float3* max);

// A routine to test whether a ray from an origin point intersects a triangle
// The routine returns a value such that if it is summed will lead to a 0 (no intersection point outside) non-zero (intersection)
float __device__ __host__ ray_triangle_intersection(float3* triangle, float3 ray_origin);

float __device__ __host__ triangel_area(float* e1, float* e2);

bool __device__ __host__ point_in_triangle(float* v0, float* v1, float* v2, float* point);

// A routine to test whether a segment from an origin point and a distance in x_axis intersects a triangle
// The routine returns a value such that if it is summed will lead to a 0 (no intersection point outside) non-zero (intersection)
float __device__ __host__ segment_triangle_intersection(float3* triangle, float3 segment_origin, float segment_length, bool *is_singular);


// A routine to test whether a segment from an origin point and a distance in x_axis intersects a triangle
// The routine returns a value such that if it is summed will lead to a 0 (no intersection point outside) non-zero (intersection)
float __device__ __host__ un_algined_segment_triangle_intersection(float3* triangle, float3 segment_origin, float3 segment_end);

// A routine to transform vertices/triangles arrays into one mesh array
void toMesh(int* triangles, float* vertices, int numtri, std::vector<float3> *mesh_vector);

// A routine to print a triangle from the triangles list given it's index
void printTriangle(int *triangles, float *vertices, int triangleIndex);

// A routine to print a triangle from a mesh given it's index
void printTriangle(float3* mesh, int triangle);

uint __device__ __host__ compute_globally_cell_state(uint *grid, uint cell_idx, float cell_intersection_value, bool is_singular);

uint __device__ compute_locally_point_state(uint origin_point_inclusion_state, float cell_intersection_value);

bool __device__ triangle_in_array(uint triangle, uint* triangles, uint triangles_length);


void __device__ get_triangles_overlap_2_cell(uint cell_id, int cell_2_id, uint* grid, uint* grid_end, uint* triangle_list, uint** triangles, uint* triangles_size);

char __device__ get_triangles_overlap_cells(uint* cell_start, uint* cell_end, uint cell_ids_size, uint* triangle_list, uint** triangles, uint* triangles_size);


void generate_random_points(float minx, float maxx, float miny, float maxy, float minz, float maxz, uint numpoints, std::vector<float3> *points_vector);

void generate_points(float minx, float maxx, float miny, float maxy, float minz, float maxz, uint numpoints, std::vector<float3> *points_vector);

void export_grid_points(uint* grid, uint gridSize, float3 min, float3 max, uint3 grid_res, float3 cell_width);

void export_test_points_as_obj_files(float3 *points, char *points_inclusion, uint numpoints);

#endif