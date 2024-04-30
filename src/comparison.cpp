#include "algo.h"
#include "utils.h"
#include "corridor_cpu.h"
#include "corridor.cuh"
#include "aabbtree_gpu.h"

#include <chrono>

#include <iostream>
#include <map>
#include <unordered_map>
#include <set>
#include <string>
#include <cassert>


// generate random variable between a and b
float generate_random_float(float a, float b) {
    assert(b > a); 
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = b - a;
    float r = random * diff;
    return a + r;
}


void load_single_organ(const std::string &organ_path, std::vector<Mymesh> &single_organ)
{
    for (fs::directory_entry& AS : fs::directory_iterator(organ_path)) 
    {
        std::string file_path = AS.path().string();
        single_organ.push_back(Mymesh(file_path));
    }

    for (auto &AS: single_organ) {
        AS.create_aabb_tree();
    }

}


float compute_intersection_volume_own_cpu(myspatial::AABBTree *p_aabbtree, Mytissue &tissue)
{
    std::vector<Point> &points = tissue.get_points();
    int count = 0;

    for (Point &p: points)
    {
        float3 p_f = make_float3(p[0], p[1], p[2]);
        if (p_aabbtree->point_inside(p_f)) count += 1;
    }
    
    return 1.0 * count / points.size();


}


int main(int argc, char **argv)
{
    //Arguments
    std::string organ_path = std::string(argv[1]);
    
    // Implementation from scratch
    Organ organ;
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

    // Using CGAL from collision detection API
    std::vector<Mymesh> organ_cgal;
    load_single_organ(organ_path, organ_cgal);


    // Fake tissue, assume AABB
    double example_d_x = 0.001;
    double example_d_y = 0.001;
    double example_d_z = 0.001;
    double tbv = example_d_x * example_d_y * example_d_z;

    int n_meshes = meshes_vector.size();

    double total_time_cgal = 0, total_time_cpu = 0;
    for (int i = 0; i < n_meshes; i++)
    {
        Mymesh &mesh_cgal = organ_cgal[i];
        Mesh &mesh = meshes_vector[i];
        MBB &bbox = mesh.bbox;

        while (true)
        {

            float c_x = generate_random_float(bbox.xmin, bbox.xmax);
            float c_y = generate_random_float(bbox.ymin, bbox.ymax);
            float c_z = generate_random_float(bbox.zmin, bbox.zmax);

            Mytissue fake_tissue(c_x, c_y, c_z, example_d_x, example_d_y, example_d_z);

            auto t1 = std::chrono::high_resolution_clock::now();
            double abs_volume = compute_intersection_volume(mesh_cgal, fake_tissue);
            double percentage1 = abs_volume / tbv;
            auto t2 = std::chrono::high_resolution_clock::now();

            auto t3 = std::chrono::high_resolution_clock::now();
            float percentage2 = compute_intersection_volume_own_cpu(p_aabbtrees[i], fake_tissue);
            auto t4 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration1 = t2 - t1;
            std::chrono::duration<double> duration2 = t4 - t3;

            total_time_cgal += duration1.count();
            total_time_cpu += duration2.count();
            std::cout << "Intersection Percentage cgal: " << percentage1 << " " << duration1.count() << "s" << ", our own: " << percentage2 << " " << duration2.count() << "s" << std::endl;
            if (percentage1 > 0)
                break; 

        }  

    }

    std::cout << "\ntotal time: -------------------" << std::endl;
    std::cout << "cgal: " << total_time_cgal << ", our own: " << total_time_cpu << std::endl;


}