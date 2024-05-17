#include "corridor_with_spatial_index.cuh"
namespace fs = boost::filesystem;

ResultContainer test_corridor_with_spatial_index_gpu(AATissue &example_tissue, std::vector<std::pair<int, float>> &result, std::vector<myspatial::AABBTree*> &p_aabbtrees, float tolerance)
{
	
	// Step 1: Copy data from host to device
	int n_collisions = result.size();

	// intermediate variable to store device memory addresses on host
	myspatial::AABBTreeCUDA *h_aabbtrees = (myspatial::AABBTreeCUDA*) malloc(n_collisions * sizeof(myspatial::AABBTreeCUDA));

	for (int i = 0; i < n_collisions; i++)
	{
		auto j = result[i].first;
		myspatial::AABBTree *p_aabbtree = p_aabbtrees[j];
		
		std::vector<myspatial::Node> &node_pool = p_aabbtree->node_pool;
		std::vector<myspatial::Triangle> &triangles = p_aabbtree->triangles_;

		myspatial::Node *h_node_pool = node_pool.data();
		myspatial::Triangle *h_triangles = triangles.data();

		int n_nodes = node_pool.size();
		int n_triangles = triangles.size();

		// Allocate memory on device and return the address to the host 
		cudaMalloc(&(h_aabbtrees[i].nodes_), n_nodes * sizeof(myspatial::Node));
		cudaMalloc(&(h_aabbtrees[i].triangles_), n_triangles * sizeof(myspatial::Triangle));

		cudaMemcpy(h_aabbtrees[i].nodes_, h_node_pool,  n_nodes * sizeof(myspatial::Node), cudaMemcpyHostToDevice);
		cudaMemcpy(h_aabbtrees[i].triangles_, h_triangles,  n_triangles * sizeof(myspatial::Triangle), cudaMemcpyHostToDevice);

		h_aabbtrees[i].n_nodes = n_nodes;
		h_aabbtrees[i].n_triangles = n_triangles;

	}
	// Copy the pointers to the device
	myspatial::AABBTreeCUDA* d_aabbtrees;
	cudaMalloc((void**)&d_aabbtrees, n_collisions * sizeof(myspatial::AABBTreeCUDA));
	cudaMemcpy(d_aabbtrees, h_aabbtrees, n_collisions * sizeof(myspatial::AABBTreeCUDA), cudaMemcpyHostToDevice);

	// Only copy intersection percentages accordingly
	std::vector<float> target_intersection_pctgs_vector;
	for (auto s: result) target_intersection_pctgs_vector.push_back(s.second);
	float *h_target_pctg = target_intersection_pctgs_vector.data();
	float *d_target_pctg;
	cudaMalloc((void**)&d_target_pctg, n_collisions * sizeof(float));
	cudaMemcpy(d_target_pctg, h_target_pctg, n_collisions * sizeof(float), cudaMemcpyHostToDevice);

	// Pass the result container
	ResultContainer result_container;
	ResultContainer *h_result_container = &result_container;
	ResultContainer *d_result_container;
	cudaMalloc((void**)&d_result_container, sizeof(ResultContainer));
	cudaMemcpy(d_result_container, h_result_container, sizeof(ResultContainer), cudaMemcpyHostToDevice);

	// Step 2: Compute candidate region
	float intersect_x_min = -1e10, intersect_y_min = -1e10, intersect_z_min = -1e10;
    float intersect_x_max = 1e10, intersect_y_max = 1e10, intersect_z_max = 1e10;  

    // size of the tissue block
    double example_d_x = example_tissue.dimension_x;
    double example_d_y = example_tissue.dimension_y;
    double example_d_z = example_tissue.dimension_z;

    double tbv = example_d_x * example_d_y * example_d_z;

    // Compute candidate region by CPU
    for (auto s: result)
    {
        // MBB &bbox = organ.meshes_vector[s.first].bbox;
		myspatial::MBB_Tri &bbox = p_aabbtrees[s.first]->get_mbb();
        intersect_x_min = std::max(intersect_x_min, bbox.xmin);
        intersect_y_min = std::max(intersect_y_min, bbox.ymin);
        intersect_z_min = std::max(intersect_z_min, bbox.zmin);

        intersect_x_max = std::min(intersect_x_max, bbox.xmax);
        intersect_y_max = std::min(intersect_y_max, bbox.ymax);
        intersect_z_max = std::min(intersect_z_max, bbox.zmax);
    }

	// Step size
    double step_x = (intersect_x_max - intersect_x_min + example_d_x) / 40.0;
    double step_y = (intersect_y_max - intersect_y_min + example_d_y) / 40.0;
    double step_z = (intersect_z_max - intersect_z_min + example_d_z) / 40.0;

    std::cout << "candidate region: " << std::endl;
    std::cout << "min x, y, z: " << intersect_x_min << " " << intersect_y_min << " " << intersect_z_min << std::endl;
    std::cout << "max x, y, z: " << intersect_x_max << " " << intersect_y_max << " " << intersect_z_max << std::endl;
    std::cout << "step size: " << step_x << " " << step_y << " " << step_z << std::endl;

	// Define Grid Configuration
	// dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
	// dim3 dimGrid(MAX_BLOCK_NUMBERS, MAX_BLOCK_NUMBERS, MAX_BLOCK_NUMBERS);
	dim3 dimBlock(512, 1, 1);
	dim3 dimGrid(64, 64, 64);

    //Launch corridor GPU
    GpuTimer timer;
	timer.Start();

	int resolution = 10;

	compute_corridor_GPU_with_spatial_index<<<dimGrid, dimBlock>>>(d_aabbtrees, d_target_pctg, n_collisions, 
                                        intersect_x_min, intersect_y_min, intersect_z_min,
                                        intersect_x_max, intersect_y_max, intersect_z_max,
                                        example_d_x, example_d_y, example_d_z,
                                        step_x, step_y, step_z, resolution, tolerance,
										d_result_container);

	// test_gpu_spatial_index<<<dimGrid, dimBlock>>>(d_aabbtrees, n_collisions);
	cudaDeviceSynchronize();
    print_if_cuda_error(__LINE__);
	timer.Stop();
	printf("\t\nKernel Time: %f msecs.\n", timer.Elapsed());

	cudaMemcpy(h_result_container, d_result_container, sizeof(ResultContainer), cudaMemcpyDeviceToHost);
	
	// why? It's so strange! why it is h_aabbtrees[i] not d_aabbtrees[i]?? still confused about cudaMalloc
	for (int i = 0; i < n_collisions; i++)
	{
		cudaFree(h_aabbtrees[i].nodes_);
		cudaFree(h_aabbtrees[i].triangles_);
	}
	
	cudaFree(d_aabbtrees);
	cudaFree(d_target_pctg);
	// cudaFree(d_result_container);

	return result_container;
	

}

__global__ void test_gpu_spatial_index(myspatial::AABBTreeCUDA* d_aabbtrees, int n_aabbtrees)
{
	int x = blockIdx.x;
    int y = blockIdx.y;
    int z = blockIdx.z;
	int i = threadIdx.x;
	int j = threadIdx.y;
	int k = threadIdx.z;

	// printf("x: %d, y: %d, z:%d\n", x, y, z);

	if (x < n_aabbtrees)
	{
		myspatial::AABBTreeCUDA &aabbtree = d_aabbtrees[x];
		myspatial::Triangle t = aabbtree.triangles_[0];
		int n_nodes = aabbtree.n_nodes;
		int n_triangles = aabbtree.n_triangles;

		// printf("x: %d, nodes: %p, triangles: %p\n", x, (void*) aabbtree.nodes_, (void*) aabbtree.triangles_);
		// printf("nodes: %d, triangles: %d, %f, %f, %f\n", n_nodes, n_triangles, t.p1.x, t.p1.y, t.p1.z);
	}
	// printf("success!\n");

}


__global__ void compute_corridor_GPU_with_spatial_index(myspatial::AABBTreeCUDA* d_aabbtrees, float* target_intersection_pctgs, uint n_meshes, 
                                        float intersect_x_min, float intersect_y_min, float intersect_z_min,
                                        float intersect_x_max, float intersect_y_max, float intersect_z_max,
                                        float example_d_x, float example_d_y, float example_d_z,
                                        float step_x, float step_y, float step_z, int resolution, float tolerance,
										ResultContainer *result_container)
{
    
    // // n_meshes is the number of collided meshes

    int x = blockIdx.x;
    int y = blockIdx.y;
    int z = blockIdx.z;

	// printf("block index: %d %d %d %d %d %d\n", x, y, z, i, j, k);
	
	if (x >= 40 || y >= 40 || z >= 40) return;

	// __shared__ float intersection_volumes[10];
	__shared__ float intersection_volumes[5][1000];
	if (threadIdx.x == 0)
	{
		// for (int i = 0; i < 10; i++) intersection_volumes[i] = 0.0;
		for (int i = 0; i < 5; i++)
			for (int j = 0; j < 1000; j++)
			intersection_volumes[i][j] = 0.0;
	}
	__syncthreads();

	float delta_x = example_d_x / resolution, delta_y = example_d_y / resolution, delta_z = example_d_z / resolution;

	// center of the tissue
	float c_x = intersect_x_min - example_d_x / 2 + step_x * x;
	float c_y = intersect_y_min - example_d_y / 2 + step_y * y;
	float c_z = intersect_z_min - example_d_z / 2 + step_z * z;

	// min of the tissue
	float min_x = c_x - example_d_x / 2;
	float min_y = c_y - example_d_y / 2;
	float min_z = c_z - example_d_z / 2;
	
	int v = 0;
	// printf("%d", count);
	
	for (int tid = threadIdx.x; tid < resolution*resolution*resolution; tid += blockDim.x)
	{
		int tmp = tid;
		int k = tmp % resolution;
		tmp /= resolution;
		int j = tmp % resolution;
		int i = tmp / resolution;

		float point_c_x = min_x + (i + 0.5) * delta_x;
		float point_c_y = min_y + (j + 0.5) * delta_y;
		float point_c_z = min_z + (k + 0.5) * delta_z;
		float3 p = make_float3(point_c_x, point_c_y, point_c_z);
				
		for (int m_idx = 0; m_idx < n_meshes; m_idx++)
		{
			// v = point_in_polyhedron(p, meshes, offset, m_idx);
			//************************attention************************
			// v = d_aabbtrees[m_idx].point_inside(p);
			v = d_aabbtrees[m_idx].point_inside(p);
			// printf("point inside: %d\n", w);
			// atomicAdd(&intersection_volumes[m_idx], v);
			intersection_volumes[m_idx][i*10*10 + j*10 + k] = v;
		}

		 
	}

	__syncthreads();

	if (threadIdx.x == 0)
	{
		bool is_in_corridor = true;
		int total_voxels = resolution * resolution * resolution;
		for (int m_idx = 0; m_idx < n_meshes; m_idx++)
		{
			float temp_intersect_volume = 0;
			for (int p_idx = 0; p_idx < 1000; p_idx++)
			{
				temp_intersect_volume += intersection_volumes[m_idx][p_idx];
			}

			// float cur_pctg = intersection_volumes[m_idx] / total_voxels;
			float cur_pctg = 1.0 * temp_intersect_volume / total_voxels;
			// printf("cur_pctg: %f\n", cur_pctg);
			if (cur_pctg > target_intersection_pctgs[m_idx] * (1 + tolerance) || cur_pctg < target_intersection_pctgs[m_idx] * (1 - tolerance))
			{
				is_in_corridor = false;
				break;
			}
		}

		int idx = x*40*40 + y*40 + z;
		if (is_in_corridor)
		{
			float3 tissue_center = make_float3(c_x, c_y, c_z);
			result_container->corridor_array[idx] = tissue_center;
			result_container->point_is_in_corridor_array[idx] = true;
		}
		else{
			result_container->point_is_in_corridor_array[idx] = false;
		}

	}
	
}


void loadAllOrganModels(std::string path, std::unordered_map<std::string, std::vector<myspatial::AABBTree*>> &total_body_gpu_with_spatial_index)
{
    for (fs::directory_entry& organ_entry: fs::directory_iterator(path)) 
    {
        Organ organ;
        std::string organ_name = organ_entry.path().stem().string();
        std::string organ_path = organ_entry.path().string();
		std::cout << organ_path << std::endl;
        loadOrganModel(organ_path, organ);

		std::vector<Mesh> &meshes_vector = organ.meshes_vector;

		// A vector of pointers of aabbtree 
		std::vector<myspatial::AABBTree*> p_aabbtrees;
		for (Mesh &mesh: meshes_vector)
		{
			std::vector<float3> &mesh_data = mesh.data;
			// Should be a vector of triangles, not vector of float3
			assert(mesh_data.size() / 3 == 0);
			std::vector<myspatial::Triangle> tri_mesh;

			for (int i = 0; i < mesh_data.size(); i+= 3)
			{
				tri_mesh.push_back(myspatial::Triangle(mesh_data[i], mesh_data[i+1], mesh_data[i+2]));

			}
			myspatial::AABBTree *aabbtree = new myspatial::AABBTree(tri_mesh);
			p_aabbtrees.push_back(aabbtree);
		}
		total_body_gpu_with_spatial_index[organ_name] = p_aabbtrees;
    }

}