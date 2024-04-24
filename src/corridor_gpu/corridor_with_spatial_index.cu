#include "corridor_with_spatial_index.cuh"


void test_corridor_with_spatial_index_gpu(AATissue &example_tissue, std::vector<std::pair<int, float>> &result, std::vector<myspatial::AABBTree*> &p_aabbtrees)
{
	// More complicated to copy an array of aabbtrees from host to GPU
	int N = result.size();
	myspatial::AABBTreeCUDA* d_aabbtrees;
	cudaMalloc((void**)&d_aabbtrees, result.size() * sizeof(myspatial::AABBTreeCUDA));

	for (auto s: result)
	{
		auto i = s.first;
		myspatial::AABBTree *p_aabbtree = p_aabbtrees[i];
		
		std::vector<myspatial::Node> &node_pool = p_aabbtree->node_pool;
		std::vector<myspatial::Triangle> &triangles = p_aabbtree->triangles_;

		myspatial::Node *h_node_pool = node_pool.data();
		myspatial::Triangle *h_triangles = triangles.data();

		int n_nodes = node_pool.size();
		int n_triangles = triangles.size();
		
		myspatial::Node *d_node;
		myspatial::Triangle *d_triangle;
		
		cudaMalloc((void**)&(d_node), n_nodes * sizeof(myspatial::Node));
		cudaMalloc((void**)&(d_triangle), n_triangles * sizeof(myspatial::Triangle));

		
		cudaMemcpy(&(d_aabbtrees[i].nodes_), &(d_node), sizeof(myspatial::Node *), cudaMemcpyHostToDevice);
		cudaMemcpy(&(d_aabbtrees[i].triangles_), &(d_triangle), sizeof(myspatial::Triangle *), cudaMemcpyHostToDevice);

		cudaMemcpy(d_aabbtrees[i].nodes_, h_node_pool, sizeof(myspatial::Node) * n_nodes, cudaMemcpyHostToDevice);
		cudaMemcpy(d_aabbtrees[i].triangles_, h_triangles, sizeof(myspatial::Triangle) * n_triangles, cudaMemcpyHostToDevice);

	}



}