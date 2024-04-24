#include "corridor_with_spatial_index.cuh"

#include <iostream>
#include <cassert>

using namespace myspatial;

int main(int argc, char **argv)
{

    // example 1, four triangles constructed by four points

    // float3 p1 = make_float3(0.0, 0.0, 0.0);
    // float3 p2 = make_float3(1.0, 0.0, 0.0);
    // float3 p3 = make_float3(0.0, 1.0, 0.0);
    // float3 p4 = make_float3(0.0, 0.0, 1.0);

    // Triangle t1(p1, p2, p3);
    // Triangle t2(p1, p2, p4);
    // Triangle t3(p1, p3, p4);
    // Triangle t4(p2, p3, p4);

    // std::vector<Triangle> triangles{t1, t2, t3, t4};

    // AABBTree *aabbtree = new AABBTree(triangles);


    // example 2
    if (argc < 2) 
    {
        std::cout << "Please provide organ path!" << std::endl;
        return 0;
    }

    std::string organ_path = std::string(argv[1]);
    Organ organ;
    loadOrganModel(organ_path, organ);
    std::vector<Mesh> &meshes_vector = organ.meshes_vector;

    // A vector of pointers of aabbtree 
    std::vector<AABBTree*> p_aabbtrees;
    for (Mesh &mesh: meshes_vector)
    {
        std::vector<float3> &mesh_data = mesh.data;
        // Should be a vector of triangles, not vector of float3
        assert(mesh_data.size() / 3 == 0);
        std::vector<Triangle> tri_mesh;

        for (int i = 0; i < mesh_data.size(); i+= 3)
        {
            tri_mesh.push_back(Triangle(mesh_data[i], mesh_data[i+1], mesh_data[i+2]));

        }
        AABBTree *aabbtree = new AABBTree(tri_mesh);
        p_aabbtrees.push_back(aabbtree);
    }

    // for (AABBTree* p_aabbtree: p_aabbtrees) std::cout << p_aabbtree->total_node << std::endl;
    
    AABBTree* p_aabbtree0 = p_aabbtrees[p_aabbtrees.size() - 1];

    // float3 p5 = make_float3(-0.011262157000601292, 0.85468000173568726, -0.083210617303848267);
    float3 p5 = make_float3(0.1, 0.1, 0.01);
    p_aabbtree0->point_inside(p5);


    std::vector<std::pair<int, float>> collision_detection_result;
    collision_detection_result.push_back(std::make_pair(3, 0.031));
    collision_detection_result.push_back(std::make_pair(4, 0.961));

    AATissue example_tissue(0, 0, 0, 0.002, 0.003, 0.002);

    test_corridor_with_spatial_index_gpu(example_tissue, collision_detection_result, p_aabbtrees);


}