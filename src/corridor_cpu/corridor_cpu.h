#pragma once

#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Real_timer.h>

#include "algo.h"
#include "utils.h"
#include <cmath>

class Mytissue: public Mymesh {

    public:
        Mytissue() = default;
        Mytissue(double c_x, double c_y, double c_z, double d_x, double d_y, double d_z); 
        
        double center_x, center_y, center_z;
        double dimension_x, dimension_y, dimension_z;
    
    public:
        std::vector<Point> &get_points();
    
    private:
        std::vector<Point> &generate_points(int resolution);
        // Surface_mesh create_mesh();

    private:
        std::vector<Point> points;


};


std::vector<Point> find_all_locations(Mymesh &my_mesh, Mytissue &example_tissue, double intersection_percentage, double tolerance);

Surface_mesh create_corridor(std::vector<Mymesh> &meshes, Mytissue &example_tissue, std::vector<double> &intersection_percnts, double tolerance);

Surface_mesh create_corridor(std::vector<Mymesh> &organ, Mytissue &example_tissue, std::vector<std::pair<int, double>> &result, double tolerance);

std::vector<Point> create_point_cloud_corridor_for_multiple_AS(std::vector<Mymesh> &meshes, Mytissue &example_tissue, std::vector<double> &intersection_percnts, double tolerance);

std::vector<Point> create_point_cloud_corridor_for_multiple_AS(std::vector<Mymesh> &organ, Mytissue &example_tissue, std::vector<std::pair<int, double>> &result, double tolerance);

Surface_mesh create_corridor_from_point_cloud(std::vector<Point> &points);

double compute_intersection_volume(Mymesh &AS, Mytissue &tissue);

double generate_pertubation(double step);