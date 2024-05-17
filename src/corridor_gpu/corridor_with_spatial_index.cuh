#pragma once
#include "corridor.cuh"
#include "aabbtree_gpu.h"


ResultContainer test_corridor_with_spatial_index_gpu(AATissue &example_tissue, std::vector<std::pair<int, float>> &result, std::vector<myspatial::AABBTree*> &p_aabbtrees, float tolerance);

__global__ void test_gpu_spatial_index(myspatial::AABBTreeCUDA* d_aabbtrees, int n_aabbtrees);

__global__ void compute_corridor_GPU_with_spatial_index(myspatial::AABBTreeCUDA* d_aabbtrees, float* target_intersection_pctgs, uint n_meshes, 
                                        float intersect_x_min, float intersect_y_min, float intersect_z_min,
                                        float intersect_x_max, float intersect_y_max, float intersect_z_max,
                                        float example_d_x, float example_d_y, float example_d_z,
                                        float step_x, float step_y, float step_z, int resolution, float tolerance,
										ResultContainer *result_container);

void loadAllOrganModels(std::string path, std::unordered_map<std::string, std::vector<myspatial::AABBTree*>> &total_body_gpu_with_spatial_index);
