#include "corridor.cuh"
namespace fs = boost::filesystem;
// namespace PMP = CGAL::Polygon_mesh_processing;

// A routin to laod a 3D object (.off) file into an array of triangle, vertices, numv (number of vertices), and numtri (number of triangles)
void offLoader(std::string path, std::vector<int> *triangles_vector, std::vector<float> *vertices_vector, int* numv, int* numtri, float3* min, float3* max) {

	std::cout << "Loading 3D object(.off) file: " << path << std::endl;
	float minx_v = 1e10, miny_v = 1e10, minz_v = 1e10;
	float maxx_v = -1e10, maxy_v = -1e10, maxz_v = -1e10;

	std::ifstream ifs(path, std::ifstream::in);
	std::string line;
	int numv_value = 0, numtri_value = 0;
	int linesReadCounter = 0;
    int cur_vertex = 0, cur_tri = 0;
    
    float x, y, z;
    int face_num_v;

	while (ifs.good() && !ifs.eof() && std::getline(ifs, line)) {
        if (linesReadCounter == 0 || line == "") continue;
        
        std::stringstream stringstream(line);
        if (linesReadCounter == 1){stringstream >> numv_value >> numtri_value; continue;}

        if (cur_vertex < numv_value)
        {
            stringstream >> x >> y >> z;
            minx_v = std::min(x, minx_v);
            maxx_v = std::max(x, maxx_v);
            miny_v = std::min(y, miny_v);
            maxy_v = std::max(y, maxy_v);
            minz_v = std::min(z, minz_v);
            maxz_v = std::max(z, maxz_v);
            vertices_vector->push_back(x);
            vertices_vector->push_back(y);
            vertices_vector->push_back(z);
            cur_vertex ++;
        }
        else if (cur_tri < numtri_value)
        {
            stringstream >> face_num_v;
            int v_index;
            for (int i = 0; i < face_num_v; i++) {
                stringstream >> v_index;
                triangles_vector->push_back(v_index);
            }
            cur_tri ++;
        }

        linesReadCounter ++;

	}// File Line Stream Loop
	printf("\n\tThe model in the file have : %d trinagles and %d vertices\n", numtri_value, numv_value);
	printf("\tThe model bounding box min. coordinate: (%f,%f,%f)\n", minx_v, miny_v, minz_v);
	printf("\tThe model bounding box max. coordinate: (%f,%f,%f)\n", maxx_v, maxy_v, maxz_v);
	printf("\n\n");

	ifs.close();
	*numv = numv_value; *numtri = numtri_value;
	min->x = minx_v; max->x = maxx_v;
	min->y = miny_v; max->y = maxy_v;
	min->z = minz_v; max->z = maxz_v;
}


// A routine to transform vertices/triangles arrays into one mesh array
void toMeshCorridor(int* triangles, float* vertices, int numtri, std::vector<float3> *mesh_vector) {

	for (int i = 0; i < numtri; i++) {

		//Add coordinates of first vertex
		int v1_index = triangles[3 * i];
		mesh_vector->push_back(make_float3(vertices[3 * v1_index], vertices[3 * v1_index + 1], vertices[3 * v1_index + 2]));

		//Add coordinates of second vertex
		int v2_index = triangles[3 * i + 1];
		mesh_vector->push_back(make_float3(vertices[3 * v2_index], vertices[3 * v2_index + 1], vertices[3 * v2_index + 2]));

		//Add coordinates of third vertex
		int v3_index = triangles[3 * i + 2];
		mesh_vector->push_back(make_float3(vertices[3 * v3_index], vertices[3 * v3_index + 1], vertices[3 * v3_index + 2]));
	}
}

void loadOrganModel(std::string path, Organ &organ) {
    
    auto &meshes_vector = organ.meshes_vector;

    int numv, numtri;
	float3 min;
	float3 max;
    uint n_meshes = 0;

    for (fs::directory_entry& AS : fs::directory_iterator(path)) 
    {
        std::string AS_file_path = AS.path().string();

        Mesh AS_mesh;
        std::vector<float3> &mesh_vector = AS_mesh.data;
	    std::vector<float> vertices_vector;
	    std::vector<int> triangles_vector;

        offLoader(path, &triangles_vector, &vertices_vector, &numv, &numtri, &min, &max);
        toMeshCorridor(triangles_vector.data(), vertices_vector.data(), numtri, &mesh_vector);
        AS_mesh.bbox = get_mbb(vertices_vector, numv);
        meshes_vector.push_back(AS_mesh);
        n_meshes ++;
    }

    organ.n_meshes = n_meshes;

}


void loadAllOrganModels(std::string path, std::unordered_map<std::string, Organ> &total_body)
{
    for (fs::directory_entry& organ_entry: fs::directory_iterator(path)) 
    {
        Organ organ;
        std::string organ_name = organ_entry.path().stem().string();
        std::string organ_path = organ_entry.path().string();
        loadOrganModel(organ_path, organ);
        total_body.insert(std::make_pair(organ_name, organ));
    }

}


//overload: create point cloud based on the collision detection result
// std::vector<Point> create_point_cloud_corridor_for_multiple_AS(Organ &organ, AATissue &example_tissue, std::vector<std::pair<int, double>> &result, double tolerance)
// {

//     std::vector<Point> center_path;
//     std::vector<Point> point_cloud;

//     double intersect_x_min = -1e10, intersect_y_min = -1e10, intersect_z_min = -1e10;
//     double intersect_x_max = 1e10, intersect_y_max = 1e10, intersect_z_max = 1e10;  

//     // size of the tissue block
//     double example_d_x = example_tissue.dimension_x;
//     double example_d_y = example_tissue.dimension_y;
//     double example_d_z = example_tissue.dimension_z;

//     double tbv = example_d_x * example_d_y * example_d_z;

//     // compute candidate region by CPU
//     for (auto s: result)
//     {
//         Mymesh &mesh = organ[s.first];
//         CGAL::Bbox_3 bbox = PMP::bbox(mesh.get_raw_mesh());
//         intersect_x_min = std::max(intersect_x_min, bbox.xmin());
//         intersect_y_min = std::max(intersect_y_min, bbox.ymin());
//         intersect_z_min = std::max(intersect_z_min, bbox.zmin());

//         intersect_x_max = std::min(intersect_x_max, bbox.xmax());
//         intersect_y_max = std::min(intersect_y_max, bbox.ymax());
//         intersect_z_max = std::min(intersect_z_max, bbox.zmax());
//     }

//     double step_x = (intersect_x_max - intersect_x_min + example_d_x) / 40.0;
//     double step_y = (intersect_y_max - intersect_y_min + example_d_y) / 40.0;
//     double step_z = (intersect_z_max - intersect_z_min + example_d_z) / 40.0;
    
//     std::cout << "candidate region: " << std::endl;
//     std::cout << "min x, y, z: " << intersect_x_min << " " << intersect_y_min << " " << intersect_z_min << std::endl;
//     std::cout << "max x, y, z: " << intersect_x_max << " " << intersect_y_max << " " << intersect_z_max << std::endl;
//     std::cout << "step size: " << step_x << " " << step_y << " " << step_z << std::endl;

//     // can be paralleled by GPU
//     for (double c_x = intersect_x_min - example_d_x / 2; c_x < intersect_x_max + example_d_x / 2; c_x += step_x)
//         for (double c_y = intersect_y_min - example_d_y / 2; c_y < intersect_y_max + example_d_y / 2; c_y += step_y)
//             for (double c_z = intersect_z_min - example_d_z / 2; c_z < intersect_z_max + example_d_z / 2; c_z += step_z)
//             {
//                 // std::cout << c_x << " " << c_y << " " << c_z << std::endl;
//                 AATissue cur_tissue(c_x, c_y, c_z, example_d_x, example_d_y, example_d_z);
                
//                 bool is_in_corridor = true;
//                 for (auto s: result)
//                 {
//                     Mymesh &mesh = organ[s.first];
//                     double example_intersection_volume = tbv * s.second;
//                     double intersection_volume = compute_intersection_volume(mesh, cur_tissue);
//                     if (std::abs(intersection_volume - example_intersection_volume) > tolerance * example_intersection_volume)
//                     {
//                         is_in_corridor = false;
//                         break;
//                     }

//                 }

//                 if (is_in_corridor)
//                 {
//                     center_path.push_back(Point(c_x, c_y, c_z));
//                     point_cloud.push_back(Point(c_x - example_d_x / 2, c_y - example_d_y / 2, c_z - example_d_z / 2));
//                     point_cloud.push_back(Point(c_x + example_d_x / 2, c_y - example_d_y / 2, c_z - example_d_z / 2));
//                     point_cloud.push_back(Point(c_x - example_d_x / 2, c_y + example_d_y / 2, c_z - example_d_z / 2));
//                     point_cloud.push_back(Point(c_x + example_d_x / 2, c_y + example_d_y / 2, c_z - example_d_z / 2));
//                     point_cloud.push_back(Point(c_x - example_d_x / 2, c_y - example_d_y / 2, c_z + example_d_z / 2));
//                     point_cloud.push_back(Point(c_x + example_d_x / 2, c_y - example_d_y / 2, c_z + example_d_z / 2));
//                     point_cloud.push_back(Point(c_x - example_d_x / 2, c_y + example_d_y / 2, c_z + example_d_z / 2));
//                     point_cloud.push_back(Point(c_x + example_d_x / 2, c_y + example_d_y / 2, c_z + example_d_z / 2));

//                 }

//             }


//     return point_cloud;

// }



__global__ void compute_corridor_GPU(float3 *meshes, uint *offset, uint n_meshes,
                                        float intersect_x_min, float intersect_y_min, float intersect_z_min,
                                        float intersect_x_max, float intersect_y_max, float intersect_z_max,
                                        float example_d_x, float example_d_y, float example_d_z,
                                        float step_x, float step_y, float step_z, int resolution)
{
    
    // n_meshes is the number of meshes
    __shared__ float intersection_volumes[500];

    int x = blockIdx.x;
    int y = blockIdx.y;
    int z = blockIdx.z;


    for (float c_x = intersect_x_min - example_d_x / 2 + step_x * x; c_x < intersect_x_max + example_d_x / 2; c_x += step_x * gridDim.x)
        for (float c_y = intersect_y_min - example_d_y / 2 + step_y * y; c_y < intersect_y_max + example_d_y / 2; c_y += step_y * gridDim.y)
            for (float c_z = intersect_z_min - example_d_z / 2 + step_z * z; c_z < intersect_z_max + example_d_z / 2; c_z += step_z * gridDim.z)
            {
                // tissue information: c_x, c_y, c_z, example_d_x, example_d_y, example_d_z;

                // Initilize intersection_volumes 
                if (threadIdx.x < n_meshes) intersection_volumes[threadIdx.x] = 0.0;
                __syncthreads();

                
                float min_x = c_x - example_d_x/2, min_y = c_y - example_d_y/2, min_z = c_z - example_d_z/2;
                float max_x = c_x + example_d_x/2, max_y = c_y + example_d_y/2, max_z = c_z + example_d_z/2; 
                float delta_x = (max_x - min_x) / resolution, delta_y = (max_y - min_y) / resolution, delta_z = (max_z - min_z) / resolution;    
                
                int thread_number_point_inside[500] = {0};

                for (int i = threadIdx.x; i < resolution; i += blockDim.x)
                    for (int j = threadIdx.y; j < resolution; j += blockDim.y)
                        for (int k = threadIdx.z; k < resolution; k += blockDim.z)
                        {
                            float point_c_x = min_x + (i + 0.5) * delta_x;
                            float point_c_y = min_y + (j + 0.5) * delta_y;
                            float point_c_z = min_z + (k + 0.5) * delta_z;
                            
                            float3 p = make_float3(point_c_x, point_c_y, point_c_z);
                            int point_result[500] = {0};
                            point_in_polyhedrons(p, meshes, offset, n_meshes, point_result);
                            
                            for (int m = 0; m < n_meshes; m++)
                                thread_number_point_inside[m] += point_result[m]; 
                        }
                
                for (int m = 0; m < n_meshes; m++)
                {
                    float thread_intersection_volume = thread_number_point_inside[m] * delta_x * delta_y * delta_z;
                    atomicAdd(&intersection_volumes[m], thread_intersection_volume);
                }

                // Wait until all threads have done their part of work
		        __syncthreads();

                // if (intersection_volumes is within a tolerance with the target)
                //     add the position to the vector;
                
                
            }

}

void __device__ __host__ point_in_polyhedrons(float3 point, float3 *meshes, uint *offset, uint n_meshes, int point_result[])
{

    for (int k = 0; k < n_meshes; k++)
    {
        float intersection_sum_per_mesh = 0;
        
        // extract the real triangle data 
        for (int i = offset[k]; i < offset[k+1]; i++)
            intersection_sum_per_mesh += ray_triangle_intersection(&meshes[3*i], point);
        
        // point inside or outside according to odd or even number of intersections
        if ((int)intersection_sum_per_mesh % 2 == 0) point_result[k] = 0;
        else point_result[k] = 1;
    }
    
}


// test corridor for all AS in each organ
void test_corridor_for_multiple_AS(AATissue &example_tissue, std::vector<std::pair<int, double>> &result, Organ &organ, double tolerance)
{

    // std::vector<Point> center_path;
    // std::vector<Point> point_cloud;

    float intersect_x_min = -1e10, intersect_y_min = -1e10, intersect_z_min = -1e10;
    float intersect_x_max = 1e10, intersect_y_max = 1e10, intersect_z_max = 1e10;  

    // size of the tissue block
    double example_d_x = example_tissue.dimension_x;
    double example_d_y = example_tissue.dimension_y;
    double example_d_z = example_tissue.dimension_z;

    double tbv = example_d_x * example_d_y * example_d_z;

    // Please start from here
    // Compute candidate region by CPU
    for (auto s: result)
    {
        MBB &bbox = organ.meshes_vector[s.first].bbox;
        intersect_x_min = std::max(intersect_x_min, bbox.xmin);
        intersect_y_min = std::max(intersect_y_min, bbox.ymin);
        intersect_z_min = std::max(intersect_z_min, bbox.zmin);

        intersect_x_max = std::min(intersect_x_max, bbox.xmax);
        intersect_y_max = std::min(intersect_y_max, bbox.ymax);
        intersect_z_max = std::min(intersect_z_max, bbox.zmax);
    }

    std::vector<uint> offsets_vector;
    std::vector<float3> collided_as;
    auto &meshes_vector = organ.meshes_vector;
    // create meshes and offsets vectors for collided meshes
    offsets_vector.push_back(0);
    // Please start from here
    for (auto s: result)
    {
        auto i = s.first;
        Mesh &mesh = meshes_vector[i];
        collided_as.insert(collided_as.end(), mesh.data.begin(), mesh.data.end());
        offsets_vector.push_back(collided_as.size());
    }

    // Number of meshes
    uint n_meshes = offsets_vector.size() - 1;
    
    // Address for the host
    float3 *h_meshes = NULL;
    uint *h_offsets = NULL;
	h_meshes = collided_as.data();
    h_offsets = offsets_vector.data();

    //Allocate memory for the device
    thrust::device_ptr<float3> d_meshes;
    thrust::device_ptr<uint> d_offsets; 
	d_meshes = thrust::device_malloc<float3>(collided_as.size());
	d_offsets = thrust::device_malloc<uint>(offsets_vector.size());
	
	// Copy from host to device
	thrust::copy(h_meshes, h_meshes + collided_as.size(), d_meshes);
	thrust::copy(h_offsets, h_offsets + offsets_vector.size(), d_offsets);

    // Step size
    double step_x = (intersect_x_max - intersect_x_min + example_d_x) / 40.0;
    double step_y = (intersect_y_max - intersect_y_min + example_d_y) / 40.0;
    double step_z = (intersect_z_max - intersect_z_min + example_d_z) / 40.0;

    std::cout << "candidate region: " << std::endl;
    std::cout << "min x, y, z: " << intersect_x_min << " " << intersect_y_min << " " << intersect_z_min << std::endl;
    std::cout << "max x, y, z: " << intersect_x_max << " " << intersect_y_max << " " << intersect_z_max << std::endl;
    std::cout << "step size: " << step_x << " " << step_y << " " << step_z << std::endl;

    //Define Grid Configuration
	dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE);
	dim3 dimGrid(MAX_BLOCK_NUMBERS, MAX_BLOCK_NUMBERS, MAX_BLOCK_NUMBERS);

    //Launch corridor GPU
    GpuTimer timer;
	timer.Start();
    int resolution = 40;

    compute_corridor_GPU<<<dimGrid, dimBlock>>>(thrust::raw_pointer_cast(d_meshes), thrust::raw_pointer_cast(d_offsets), n_meshes,
                                        intersect_x_min, intersect_y_min, intersect_z_min,
                                        intersect_x_max, intersect_y_max, intersect_z_max,
                                        example_d_x, example_d_y, example_d_z,
                                        step_x, step_y, step_z, resolution);
    
    cudaDeviceSynchronize();
    print_if_cuda_error(__LINE__);
	timer.Stop();
	printf("\t\n Kernel Time: %f msecs.\n", timer.Elapsed());


    // copy from device to host
	// should be added in the future

	cudaFree(thrust::raw_pointer_cast(d_meshes));
	cudaFree(thrust::raw_pointer_cast(d_offsets));

}


MBB get_mbb(std::vector<float>& vertices_vector, int numv)
{
    float x, y, z;
    float minx_v = 1e10, miny_v = 1e10, minz_v = 1e10;
	float maxx_v = -1e10, maxy_v = -1e10, maxz_v = -1e10;

    for (int i = 0; i < numv; i++)
    {
        x = vertices_vector[3*i];
        y = vertices_vector[3*i + 1];
        z = vertices_vector[3*i + 2];

        minx_v = std::min(x, minx_v);
        maxx_v = std::max(x, maxx_v);
        miny_v = std::min(y, miny_v);
        maxy_v = std::max(y, maxy_v);
        minz_v = std::min(z, minz_v);
        maxz_v = std::max(z, maxz_v);
    }

    return MBB(minx_v, miny_v, minz_v, maxx_v, maxy_v, maxz_v);

}