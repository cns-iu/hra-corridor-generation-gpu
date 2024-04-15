#include "aabbtree.h"
#include "corridor.cuh"

#include <iostream>
#include <cassert>

int main(int argc, char **argv)
{

    // example 1, four triangles constructed by four points

    float3 p1 = make_float3(0.0, 0.0, 0.0);
    float3 p2 = make_float3(1.0, 0.0, 0.0);
    float3 p3 = make_float3(0.0, 1.0, 0.0);
    float3 p4 = make_float3(0.0, 0.0, 1.0);

    Triangle t1(p1, p2, p3);
    Triangle t2(p1, p2, p4);
    Triangle t3(p1, p3, p4);
    Triangle t4(p2, p3, p4);

    std::vector<Triangle> triangles{t1, t2, t3, t4};

    AABBTree *aabbtree = new AABBTree(triangles);


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
    }
    

}