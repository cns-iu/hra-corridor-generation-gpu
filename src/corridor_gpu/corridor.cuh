#pragma once

//Misc Libs
#include <stdlib.h>
#include <stdio.h>
#include<math.h>
#include <string>
#include<time.h>

// Point In Polynomial Helpers
// #include <pip_helpers.cuh>

//Cuda Libs
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
// #include "device_functions.h"
#include "GpuTimer.h"

//Thrust
#include <thrust/scan.h>
#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/fill.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <unordered_map>
#include <boost/filesystem.hpp>

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


#define BLOCK_SIZE 1024
#define NUM_POINTS 100000
#define MAX_BLOCK_NUMBERS 65535


// Axis-Aligned Tissue
class AATissue{
    public:
        AATissue() = default;
        AATissue(float c_x, float c_y, float c_z, float d_x, float d_y, float d_z):
        center_x(c_x), center_y(c_y), center_z(c_z), dimension_x(d_x), dimension_y(d_y), dimension_z(d_z) {}

    public:
        float center_x, center_y, center_z;
        float dimension_x, dimension_y, dimension_z;

};

class MBB{
    public:
        MBB() = default;
        MBB(float x_min, float y_min, float z_min, float x_max, float y_max, float z_max):
        xmin(x_min), ymin(y_min), zmin(z_min), xmax(x_max), ymax(y_max), zmax(z_max) {}

    public:
        float xmin, ymin, zmin, xmax, ymax, zmax;

};

class Mesh{
    public:
        std::vector<float3> data;
        std::string label;
        MBB bbox;
    
};

class Organ{
    public:
        // meshes_vector is the concatenation of a 3d points. Three 3d points form a face. 
        // (offset[i+1] - offset[i]) * 3 form the i-th mesh
        std::vector<Mesh> meshes_vector;
        // std::unordered_map<Mesh> meshes_map;
        uint n_meshes;
};

class ResultContainer{
    public:
        float3 corridor_array[40*40*40];
        bool point_is_in_corridor_array[40*40*40];
        int corridor_point_idx = 0;
};

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

// Load 3D object in OFF
void offLoader(const char *path, std::vector<int> *triangles_vector, std::vector<float> *vertices_vector, int* numv, int* numtri, float3* min, float3* max);

// Convert triangle vector and vertice vector into one mesh vector
void toMeshCorridor(int* triangles, float* vertices, int numtri, std::vector<float3> *mesh_vector);

// Load multiple meshes into one organ
void loadOrganModel(std::string path, Organ &organ);

// Load all the organs
void loadAllOrganModels(std::string path, std::unordered_map<std::string, Organ> &total_body);
int add(int a, int b);


__global__ void compute_corridor_GPU(float3 *meshes, uint *offset, float* target_intersection_pctgs, uint n_meshes, 
                                        float intersect_x_min, float intersect_y_min, float intersect_z_min,
                                        float intersect_x_max, float intersect_y_max, float intersect_z_max,
                                        float example_d_x, float example_d_y, float example_d_z,
                                        float step_x, float step_y, float step_z, int resolution, float tolerance,
                                        ResultContainer *result_container);

void __device__ __host__ point_in_polyhedrons(float3 point, float3 *meshes, uint *offset, uint n_meshes, int point_result[]);

int __device__ __host__ point_in_polyhedron(float3 point, float3 *meshes, uint *offset, int m_index);

// void test_corridor_for_multiple_AS(std::vector<Mymesh> &organ, AATissue &example_tissue, std::vector<std::pair<int, double>> &result, double tolerance);

ResultContainer test_corridor_for_multiple_AS(AATissue &example_tissue, std::vector<std::pair<int, float>> &result, Organ &organ, float tolerance);

MBB get_mbb(std::vector<float>& vertices_vector, int numv);


