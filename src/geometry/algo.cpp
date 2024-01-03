#include "algo.h"

std::vector<std::pair<std::string, double>> collision_detection_single_tissue(std::vector<Mymesh> &organ, Mymesh &tissue, std::vector<Point> &points)
{
    auto aabbtree_tissue = tissue.get_aabb_tree();
    std::vector<std::pair<std::string, double>> result;

    // std::cout << "organ size: " << organ.size() << std::endl;
    for (int i = 0; i < organ.size(); i++)
    {
        auto AS = organ[i];
        
        if (!AS.is_surface) {std::cout << "not surface: " << AS.label << std::endl; continue;}

        auto aabbtree_AS = AS.get_aabb_tree();
        if (aabbtree_AS->do_intersect(*aabbtree_tissue))
        {
            
            if (AS.is_closed) 
            {
                double percentage = AS.percentage_points_inside(points);
                result.push_back({AS.label, percentage});
            }
            else
                result.push_back({AS.label, -1.0});
        }
        else
        {

            Surface_mesh &mesh = tissue.get_raw_mesh();
            
            // the tissue block is wholely inside the anatomical structure. 
            bool is_contain_1 = true;
            for (auto vd: mesh.vertices())
            {
                Point p = mesh.point(vd);
                if (!AS.point_inside(p)) 
                {
                    is_contain_1 = false; 
                    break;
                }
            }

            // the anatomical structure is wholely inside the tissue block, still use the voxel-based algorithm, can be simplified to use the volume of the anatomical structure. 
            bool is_contain_2 = true;
            Surface_mesh &AS_raw_mesh = AS.get_raw_mesh();

            for (auto vd: AS_raw_mesh.vertices())
            {
                Point p = AS_raw_mesh.point(vd);
                
                if (!tissue.point_inside(p))
                    is_contain_2 = false;
                break;
            }

    
            if (is_contain_1) 
                result.push_back({AS.label, 1.0});
            else if (is_contain_2)
            {
                double percentage = AS.percentage_points_inside(points);
                result.push_back({AS.label, percentage});
            }
        }
    }

    return result;
}

std::vector<std::string> collision_detection_bb(std::vector<Mymesh> &organ, Mymesh &tissue)
{
    std::vector<std::string> result;

    CGAL::Bbox_3 bb_tissue = PMP::bbox(tissue.get_raw_mesh());

    for (auto &AS: organ)
    {
        CGAL::Bbox_3 bb_AS = PMP::bbox(AS.get_raw_mesh());
        if (CGAL::do_overlap(bb_AS, bb_tissue))
            result.push_back(AS.label);
    }

    return result;
}

std::vector<std::pair<int, double>> collision_detection_single_tissue_2(std::vector<Mymesh> &organ, Mymesh &tissue, std::vector<Point> &points)
{
    auto aabbtree_tissue = tissue.get_aabb_tree();
    std::vector<std::pair<int, double>> result;

    // std::cout << "organ size: " << organ.size() << std::endl;
    for (int i = 0; i < organ.size(); i++)
    {
        auto &AS = organ[i];
        
        if (!AS.is_surface) {std::cout << "not surface: " << AS.label << std::endl; continue;}

        auto aabbtree_AS = AS.get_aabb_tree();
        if (aabbtree_AS->do_intersect(*aabbtree_tissue))
        {
            
            if (AS.is_closed) 
            {
                double percentage = AS.percentage_points_inside(points);
                result.push_back({i, percentage});
            }
            else
                result.push_back({i, -1.0});
        }
        else
        {

            Surface_mesh &mesh = tissue.get_raw_mesh();
            
            // the tissue block is wholely inside the anatomical structure. 
            bool is_contain_1 = true;
            for (auto vd: mesh.vertices())
            {
                Point p = mesh.point(vd);
                if (!AS.point_inside(p)) 
                {
                    is_contain_1 = false; 
                    break;
                }
            }

            // the anatomical structure is wholely inside the tissue block, still use the voxel-based algorithm, can be simplified to use the volume of the anatomical structure. 
            bool is_contain_2 = true;
            Surface_mesh &AS_raw_mesh = AS.get_raw_mesh();

            for (auto vd: AS_raw_mesh.vertices())
            {
                Point p = AS_raw_mesh.point(vd);
                
                if (!tissue.point_inside(p))
                    is_contain_2 = false;
                break;
            }

    
            if (is_contain_1) 
                result.push_back({i, 1.0});
            else if (is_contain_2)
            {
                double percentage = AS.percentage_points_inside(points);
                result.push_back({i, percentage});
            }
        }
    }

    return result;
}
