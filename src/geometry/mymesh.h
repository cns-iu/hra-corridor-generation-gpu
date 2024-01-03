#ifndef MYMESH_H_
#define MYMESH_H_

// #include<CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/Surface_mesh.h>
#include<CGAL/polygon_mesh_processing.h>
#include<CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include<fstream>
#include<iostream>
#include<vector>
#include<boost/filesystem.hpp>
#include "omp.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Triangle_3 Triangle;

typedef CGAL::Surface_mesh<Point> Surface_mesh;

typedef CGAL::AABB_face_graph_triangle_primitive<Surface_mesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Surface_mesh::Face_range Face_range;
typedef Surface_mesh::Halfedge_index halfedge_descriptor;
typedef Surface_mesh::Vertex_index vertex_descriptor;
typedef Surface_mesh::Face_index face_descriptor;
typedef CGAL::Side_of_triangle_mesh<Surface_mesh, Kernel> Point_inside;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace fs = boost::filesystem;

class Mymesh{

    public:
        Mymesh(){};
        Mymesh(const std::string &file_path);
        Mymesh(const Surface_mesh sm);
        Mymesh(double c_x, double c_y, double c_z, double d_x, double d_y, double d_z);

        ~Mymesh(){};

        bool load_from_off(const std::string &file_path);
        bool triangulate_mesh();
        bool point_inside(Point &query);
        double percentage_points_inside(std::vector<Point> &query);

        void create_aabb_tree();
        Surface_mesh& get_raw_mesh();
        Tree* get_aabb_tree();

    private:
        Surface_mesh mesh;

    public:
        std::shared_ptr<Tree> aabbTree;
        std::string label;
        bool is_closed;
        bool is_surface;
        double volume = 0.0;
        std::string to_wkt();
};

#endif