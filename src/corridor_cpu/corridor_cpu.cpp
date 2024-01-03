#include "corridor_cpu.h"


Mytissue::Mytissue(double c_x, double c_y, double c_z, double d_x, double d_y, double d_z)
:Mymesh(c_x, c_y, c_z, d_x, d_y, d_z),
center_x(c_x), center_y(c_y), center_z(c_z), dimension_x(d_x), dimension_y(d_y), dimension_z(d_z)
{

    create_aabb_tree();
    generate_points(10);

}


std::vector<Point> &Mytissue::get_points()
{
    return this->points;
}


// Surface_mesh Mytissue::create_mesh()
// {

//     double min_x = center_x - dimension_x/2, min_y = center_y - dimension_y/2, min_z = center_z - dimension_z/2;
//     double max_x = center_x + dimension_x/2, max_y = center_y + dimension_y/2, max_z = center_z + dimension_z/2;

//     Point v000(min_x, min_y, min_z);
//     Point v100(max_x, min_y, min_z);
//     Point v010(min_x, max_y, min_z);
//     Point v001(min_x, min_y, max_z);
//     Point v110(max_x, max_y, min_z);
//     Point v101(max_x, min_y, max_z);
//     Point v011(min_x, max_y, max_z);
//     Point v111(max_x, max_y, max_z);

//     std::vector<Point> vertices = {v000, v100, v110, v010, v001, v101, v111, v011};
//     std::vector<vertex_descriptor> vd;

//     Surface_mesh tissue_mesh;
//     for (auto &p: vertices)
//     {
//         vertex_descriptor u = tissue_mesh.add_vertex(p);
//         vd.push_back(u);
//     } 

//     tissue_mesh.add_face(vd[3], vd[2], vd[1], vd[0]);
//     tissue_mesh.add_face(vd[4], vd[5], vd[6], vd[7]);
//     tissue_mesh.add_face(vd[4], vd[7], vd[3], vd[0]);
//     tissue_mesh.add_face(vd[1], vd[2], vd[6], vd[5]);
//     tissue_mesh.add_face(vd[0], vd[1], vd[5], vd[4]);
//     tissue_mesh.add_face(vd[2], vd[3], vd[7], vd[6]);

//     return tissue_mesh;


// }

std::vector<Point> &Mytissue::generate_points(int resolution=10)
{

    double min_x = center_x - dimension_x/2, min_y = center_y - dimension_y/2, min_z = center_z - dimension_z/2;
    double max_x = center_x + dimension_x/2, max_y = center_y + dimension_y/2, max_z = center_z + dimension_z/2; 
    double delta_x = (max_x - min_x) / resolution, delta_y = (max_y - min_y) / resolution, delta_z = (max_z - min_z) / resolution;    
  

    for (int i = 0; i < resolution; i++)
        for (int j = 0; j < resolution; j++)
            for (int k = 0; k < resolution; k++)
            {
                double c_x = min_x + (i + 0.5) * delta_x;
                double c_y = min_y + (j + 0.5) * delta_y;
                double c_z = min_z + (k + 0.5) * delta_z;
                
                Point p(c_x, c_y, c_z);
                points.push_back(p); 
            }
    
    return points;
}


std::vector<Point> find_all_locations(Mymesh &my_mesh, Mytissue &example_tissue, double intersection_percentage, double tolerance)
{

    std::vector<Point> center_path;

    // double example_intersection_volume = compute_intersection_volume(my_mesh, example_tissue);
    double example_intersection_volume = example_tissue.dimension_x * example_tissue.dimension_y * example_tissue.dimension_z * intersection_percentage;

    double step_x = example_tissue.dimension_x / 5;
    double step_y = example_tissue.dimension_y / 5;
    double step_z = example_tissue.dimension_z / 5;

    double example_d_x = example_tissue.dimension_x;
    double example_d_y = example_tissue.dimension_y;
    double example_d_z = example_tissue.dimension_z;

    Surface_mesh mesh = my_mesh.get_raw_mesh();
    CGAL::Bbox_3 bbox = PMP::bbox(mesh);

    double x_min = bbox.xmin();
    double y_min = bbox.ymin();
    double z_min = bbox.zmin();
    double x_max = bbox.xmax();
    double y_max = bbox.ymax();
    double z_max = bbox.zmax();

    for (double c_x = x_min - example_d_x / 2; c_x < x_max + example_d_x / 2; c_x += step_x)
        for (double c_y = y_min - example_d_y / 2; c_y < y_max + example_d_y / 2; c_y += step_y)
            for (double c_z = z_min - example_d_z / 2; c_z < z_max + example_d_z / 2; c_z += step_z)
            {
                Mytissue cur_tissue(c_x, c_y, c_z, example_d_x, example_d_y, example_d_z);
                double intersection_volume = compute_intersection_volume(my_mesh, cur_tissue);
                
                if (std::abs(intersection_volume - example_intersection_volume) <= tolerance * example_intersection_volume)
                    center_path.push_back(Point(c_x, c_y, c_z));

            }
    
    return center_path;
    

}

Surface_mesh create_corridor(std::vector<Mymesh> &meshes, Mytissue &example_tissue, std::vector<double> &intersection_percnts, double tolerance)
{

    std::vector<Point> point_cloud = create_point_cloud_corridor_for_multiple_AS(meshes, example_tissue, intersection_percnts, tolerance);
    return create_corridor_from_point_cloud(point_cloud);

}

// overload: create corridor based on the collision detection result
Surface_mesh create_corridor(std::vector<Mymesh> &organ, Mytissue &example_tissue, std::vector<std::pair<int, double>> &result, double tolerance)
{
    std::vector<Point> point_cloud = create_point_cloud_corridor_for_multiple_AS(organ, example_tissue, result, tolerance);
    return create_corridor_from_point_cloud(point_cloud);

}

double generate_pertubation(double step)
{
    return step * ((double)rand()/RAND_MAX * 2.0 - 1.0);
}

std::vector<Point> create_point_cloud_corridor_for_multiple_AS(std::vector<Mymesh> &meshes, Mytissue &example_tissue, std::vector<double> &intersection_percnts, double tolerance)
{

    std::vector<Point> center_path;
    std::vector<Point> point_cloud;

    double intersect_x_min = -1e10, intersect_y_min = -1e10, intersect_z_min = -1e10;
    double intersect_x_max = 1e10, intersect_y_max = 1e10, intersect_z_max = 1e10;  

    // double step_x = example_tissue.dimension_x / 4.9;
    // double step_y = example_tissue.dimension_y / 4.9;
    // double step_z = example_tissue.dimension_z / 4.9;

    double example_d_x = example_tissue.dimension_x;
    double example_d_y = example_tissue.dimension_y;
    double example_d_z = example_tissue.dimension_z;

    double tbv = example_d_x * example_d_y * example_d_z;

    for (Mymesh &mesh: meshes)
    {
        CGAL::Bbox_3 bbox = PMP::bbox(mesh.get_raw_mesh());
        intersect_x_min = std::max(intersect_x_min, bbox.xmin());
        intersect_y_min = std::max(intersect_y_min, bbox.ymin());
        intersect_z_min = std::max(intersect_z_min, bbox.zmin());

        intersect_x_max = std::min(intersect_x_max, bbox.xmax());
        intersect_y_max = std::min(intersect_y_max, bbox.ymax());
        intersect_z_max = std::min(intersect_z_max, bbox.zmax());

    }

    double step_x = (intersect_x_max - intersect_x_min + example_d_x) / 20.0;
    double step_y = (intersect_y_max - intersect_y_min + example_d_y) / 20.0;
    double step_z = (intersect_z_max - intersect_z_min + example_d_z) / 20.0;

    std::cout << "min x, y, z: " << intersect_x_min << " " << intersect_y_min << " " << intersect_z_min << std::endl;
    std::cout << "max x, y, z: " << intersect_x_max << " " << intersect_y_max << " " << intersect_z_max << std::endl;
    std::cout << "step size: " << step_x << " " << step_y << " " << step_z << std::endl;
    
    for (double c_x = intersect_x_min - example_d_x / 2; c_x < intersect_x_max + example_d_x / 2; c_x += step_x)
        for (double c_y = intersect_y_min - example_d_y / 2; c_y < intersect_y_max + example_d_y / 2; c_y += step_y)
            for (double c_z = intersect_z_min - example_d_z / 2; c_z < intersect_z_max + example_d_z / 2; c_z += step_z)
            {
                // std::cout << c_x << " " << c_y << " " << c_z << std::endl;
                Mytissue cur_tissue(c_x, c_y, c_z, example_d_x, example_d_y, example_d_z);
                
                bool is_in_corridor = true;
                for (int i = 0; i < meshes.size(); i++)
                {
                    Mymesh &mesh = meshes[i];
                    double example_intersection_volume = tbv * intersection_percnts[i];

                    double intersection_volume = compute_intersection_volume(mesh, cur_tissue);
                    // std::cout << i << " example intersection volume: " << example_intersection_volume << " " << intersection_volume << std::endl;
                    if (std::abs(intersection_volume - example_intersection_volume) > tolerance * example_intersection_volume)
                    {
                        is_in_corridor = false;
                        break;
                    }

                }

                if (is_in_corridor)
                {
                    center_path.push_back(Point(c_x, c_y, c_z));
                    std::cout << generate_pertubation(step_x) << " " << generate_pertubation(step_y) << " " << generate_pertubation(step_z) << std::endl;
                    point_cloud.push_back(Point(c_x - example_d_x / 2 + generate_pertubation(step_x), c_y - example_d_y / 2 + generate_pertubation(step_y), c_z - example_d_z / 2 + generate_pertubation(step_z)));
                    point_cloud.push_back(Point(c_x + example_d_x / 2 + generate_pertubation(step_x), c_y - example_d_y / 2 + generate_pertubation(step_y), c_z - example_d_z / 2 + generate_pertubation(step_z)));
                    point_cloud.push_back(Point(c_x - example_d_x / 2 + generate_pertubation(step_x), c_y + example_d_y / 2 + generate_pertubation(step_y), c_z - example_d_z / 2 + generate_pertubation(step_z)));
                    point_cloud.push_back(Point(c_x + example_d_x / 2 + generate_pertubation(step_x), c_y + example_d_y / 2 + generate_pertubation(step_y), c_z - example_d_z / 2 + generate_pertubation(step_z)));
                    point_cloud.push_back(Point(c_x - example_d_x / 2 + generate_pertubation(step_x), c_y - example_d_y / 2 + generate_pertubation(step_y), c_z + example_d_z / 2 + generate_pertubation(step_z)));
                    point_cloud.push_back(Point(c_x + example_d_x / 2 + generate_pertubation(step_x), c_y - example_d_y / 2 + generate_pertubation(step_y), c_z + example_d_z / 2 + generate_pertubation(step_z)));
                    point_cloud.push_back(Point(c_x - example_d_x / 2 + generate_pertubation(step_x), c_y + example_d_y / 2 + generate_pertubation(step_y), c_z + example_d_z / 2 + generate_pertubation(step_z)));
                    point_cloud.push_back(Point(c_x + example_d_x / 2 + generate_pertubation(step_x), c_y + example_d_y / 2 + generate_pertubation(step_y), c_z + example_d_z / 2 + generate_pertubation(step_z)));

                }

            }


    return point_cloud;

}

//overload: create point cloud based on the collision detection result
std::vector<Point> create_point_cloud_corridor_for_multiple_AS(std::vector<Mymesh> &organ, Mytissue &example_tissue, std::vector<std::pair<int, double>> &result, double tolerance)
{

    std::vector<Point> center_path;
    std::vector<Point> point_cloud;

    double intersect_x_min = -1e10, intersect_y_min = -1e10, intersect_z_min = -1e10;
    double intersect_x_max = 1e10, intersect_y_max = 1e10, intersect_z_max = 1e10;  

    // double step_x = example_tissue.dimension_x / 4.9;
    // double step_y = example_tissue.dimension_y / 4.9;
    // double step_z = example_tissue.dimension_z / 4.9;

    double example_d_x = example_tissue.dimension_x;
    double example_d_y = example_tissue.dimension_y;
    double example_d_z = example_tissue.dimension_z;

    double tbv = example_d_x * example_d_y * example_d_z;

    for (auto s: result)
    {
        Mymesh &mesh = organ[s.first];
        CGAL::Bbox_3 bbox = PMP::bbox(mesh.get_raw_mesh());
        intersect_x_min = std::max(intersect_x_min, bbox.xmin());
        intersect_y_min = std::max(intersect_y_min, bbox.ymin());
        intersect_z_min = std::max(intersect_z_min, bbox.zmin());

        intersect_x_max = std::min(intersect_x_max, bbox.xmax());
        intersect_y_max = std::min(intersect_y_max, bbox.ymax());
        intersect_z_max = std::min(intersect_z_max, bbox.zmax());
    }

    double step_x = (intersect_x_max - intersect_x_min + example_d_x) / 40.0;
    double step_y = (intersect_y_max - intersect_y_min + example_d_y) / 40.0;
    double step_z = (intersect_z_max - intersect_z_min + example_d_z) / 40.0;
    
    std::cout << "min x, y, z: " << intersect_x_min << " " << intersect_y_min << " " << intersect_z_min << std::endl;
    std::cout << "max x, y, z: " << intersect_x_max << " " << intersect_y_max << " " << intersect_z_max << std::endl;
    std::cout << "step size: " << step_x << " " << step_y << " " << step_z << std::endl;

    
    for (double c_x = intersect_x_min - example_d_x / 2; c_x < intersect_x_max + example_d_x / 2; c_x += step_x)
        for (double c_y = intersect_y_min - example_d_y / 2; c_y < intersect_y_max + example_d_y / 2; c_y += step_y)
            for (double c_z = intersect_z_min - example_d_z / 2; c_z < intersect_z_max + example_d_z / 2; c_z += step_z)
            {
                // std::cout << c_x << " " << c_y << " " << c_z << std::endl;
                Mytissue cur_tissue(c_x, c_y, c_z, example_d_x, example_d_y, example_d_z);
                
                bool is_in_corridor = true;
                for (auto s: result)
                {
                    Mymesh &mesh = organ[s.first];
                    double example_intersection_volume = tbv * s.second;
                    double intersection_volume = compute_intersection_volume(mesh, cur_tissue);
                    if (std::abs(intersection_volume - example_intersection_volume) > tolerance * example_intersection_volume)
                    {
                        is_in_corridor = false;
                        break;
                    }

                }

                if (is_in_corridor)
                {
                    center_path.push_back(Point(c_x, c_y, c_z));
                    point_cloud.push_back(Point(c_x - example_d_x / 2, c_y - example_d_y / 2, c_z - example_d_z / 2));
                    point_cloud.push_back(Point(c_x + example_d_x / 2, c_y - example_d_y / 2, c_z - example_d_z / 2));
                    point_cloud.push_back(Point(c_x - example_d_x / 2, c_y + example_d_y / 2, c_z - example_d_z / 2));
                    point_cloud.push_back(Point(c_x + example_d_x / 2, c_y + example_d_y / 2, c_z - example_d_z / 2));
                    point_cloud.push_back(Point(c_x - example_d_x / 2, c_y - example_d_y / 2, c_z + example_d_z / 2));
                    point_cloud.push_back(Point(c_x + example_d_x / 2, c_y - example_d_y / 2, c_z + example_d_z / 2));
                    point_cloud.push_back(Point(c_x - example_d_x / 2, c_y + example_d_y / 2, c_z + example_d_z / 2));
                    point_cloud.push_back(Point(c_x + example_d_x / 2, c_y + example_d_y / 2, c_z + example_d_z / 2));

                }

            }


    return point_cloud;

}


Surface_mesh create_corridor_from_point_cloud(std::vector<Point> &points)
{
    Surface_mesh wrap;
    
    std::cout << "points size: " << points.size() << std::endl;
    if (points.size() == 0)
        return wrap;

    // Compute the alpha and offset values
    const double relative_alpha = 10.;
    const double relative_offset = 300.;
    
    CGAL::Bbox_3 bbox = CGAL::bbox_3(std::cbegin(points), std::cend(points));
    const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                        CGAL::square(bbox.ymax() - bbox.ymin()) +
                                        CGAL::square(bbox.zmax() - bbox.zmin()));

    const double alpha = diag_length / relative_alpha;
    const double offset = diag_length / relative_offset;
    std::cout << "absolute alpha = " << alpha << " absolute offset = " << offset << std::endl;
    // Construct the wrap
    CGAL::Real_timer t;
    
    t.start();
    std::cout << "size of point cloud: " << points.size() << std::endl;
    CGAL::alpha_wrap_3(points, alpha, offset, wrap);
    t.stop();

    return wrap;

}


// compute intersection volume between a mesh and an axis-aligned tissue 
double compute_intersection_volume(Mymesh &AS, Mytissue &tissue)
{

    auto aabbtree_AS = AS.get_aabb_tree();
    auto aabbtree_tissue = tissue.get_aabb_tree();

    double percentage = 0.0;
    if (aabbtree_AS->do_intersect(*aabbtree_tissue))
    {
        percentage = AS.percentage_points_inside(tissue.get_points());
    }
    else
    {
        Surface_mesh &tissue_raw_mesh = tissue.get_raw_mesh();
        
        // the tissue block is wholely inside the anatomical structure. 
        bool is_contain_1 = true;
        for (auto vd: tissue_raw_mesh.vertices())
        {
            Point p = tissue_raw_mesh.point(vd);

            if (!AS.point_inside(p)) is_contain_1 = false;
            break;

        }


        // the anatomical structure is wholely inside the tissue block, still use the voxel-based algorithm, can be simplified to use the volume of the anatomical structure. 
        bool is_contain_2 = true;
        Surface_mesh &AS_raw_mesh = AS.get_raw_mesh();

        for (auto vd: AS_raw_mesh.vertices())
        {
            Point p = AS_raw_mesh.point(vd);

            if (!tissue.point_inside(p)) is_contain_2 = false;
            break;
        
        }

        if (is_contain_1) percentage = 1.0;
        else if (is_contain_2) percentage = AS.percentage_points_inside(tissue.get_points());
    }

    double volume = percentage * tissue.dimension_x * tissue.dimension_y * tissue.dimension_z;
    
    return volume;
    
}