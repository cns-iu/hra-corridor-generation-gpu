#include "corridor.cuh"
#include "aabbtree_gpu.h"


void test_corridor_with_spatial_index_gpu(AATissue &example_tissue, std::vector<std::pair<int, float>> &result, std::vector<myspatial::AABBTree*> &p_aabbtrees);


