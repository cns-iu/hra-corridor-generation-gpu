#pragma once

//Misc Libs
#include <stdlib.h>
#include <stdio.h>
#include<math.h>
#include <string>
#include<time.h>

// Point In Polynomial Helpers
#include "pip_helpers.cuh"

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

//CGAL
// #include "mymesh.h"


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


// Load 3D object in OFF
void offLoader(const char *path, std::vector<int> *triangles_vector, std::vector<float> *vertices_vector, int* numv, int* numtri, float3* min, float3* max);

// Convert triangle vector and vertice vector into one mesh vector
void toMeshCorridor(int* triangles, float* vertices, int numtri, std::vector<float3> *mesh_vector);

// Load multiple meshes into one organ
void loadOrganModel(std::string path, Organ &organ);

// Load all the organs
void loadAllOrganModels(std::string path, std::unordered_map<std::string, Organ> &total_body);


__global__ void compute_corridor_GPU(float3 *meshes, uint *offset, uint n_meshes, float3 *point_cloud, 
                                        float intersect_x_min, float intersect_y_min, float intersect_z_min,
                                        float intersect_x_max, float intersect_y_max, float intersect_z_max,
                                        float example_d_x, float example_d_y, float example_d_z,
                                        float step_x, float step_y, float step_z, int resolution);

void __device__ __host__ point_in_polyhedrons(float3 point, float3 *meshes, uint *offset, uint n_meshes, int point_result[]);


// void test_corridor_for_multiple_AS(std::vector<Mymesh> &organ, AATissue &example_tissue, std::vector<std::pair<int, double>> &result, double tolerance);

void test_corridor_for_multiple_AS(AATissue &example_tissue, std::vector<std::pair<int, double>> &result, double tolerance);

MBB get_mbb(std::vector<float>& vertices_vector, int numv);


