#include "mymesh.h"

namespace PMP = CGAL::Polygon_mesh_processing;
namespace fs = boost::filesystem;


bool Mymesh::load_from_off(const std::string &file_path) {
    
    std::ifstream input(file_path);

    std::string AS_name = fs::path(file_path).stem().string();

    if (!input ||!(input >> this->mesh) || this->mesh.is_empty()){
        std::cerr << file_path <<  " Not a valid mesh" << std::endl;
        this->is_surface = false;
        this->is_closed = false;
        return false;
    }

    this->label = AS_name;
    this->is_surface = true;
    if (CGAL::is_closed(this->mesh)){
        this->is_closed = true;
        this->volume = PMP::volume(this->mesh) * 1e9;
    }
    return true;

}

bool Mymesh::triangulate_mesh() {

    PMP::triangulate_faces(this->mesh);

    //confirm all faces are triangles

    for (boost::graph_traits<Surface_mesh>::face_descriptor fit : faces(mesh)) {
        
        if (next(next(halfedge(fit, mesh), mesh), mesh) 
                != prev(halfedge(fit, mesh), mesh)){
            std::cerr << "Error: non-triangular face left in mesh." << std::endl;

            return false;
        }
    }

    return true;

}

bool Mymesh::point_inside(Point& query)
{
    Point_inside inside_tester(*aabbTree);
    return inside_tester(query) == CGAL::ON_BOUNDED_SIDE;
}

double Mymesh::percentage_points_inside(std::vector<Point> &query)
{
    int cnt = 0;
    Point_inside inside_tester(*aabbTree);

    // #pragma omp parallel for reduction(+:cnt)
    for (auto &point: query)
        if (inside_tester(point) == CGAL::ON_BOUNDED_SIDE) cnt++;
    
    double percentage = 1.0 * cnt / query.size();

    return percentage;

}

void Mymesh::create_aabb_tree() {

        // std::unique_ptr<Tree> tree = std::make_unique<Tree> (faces(mesh).first, faces(mesh).second, mesh);
        std::shared_ptr<Tree> tree = std::make_shared<Tree> (faces(mesh).first, faces(mesh).second, mesh);
        // Tree* tree = new Tree(faces(mesh).first, faces(mesh).second, mesh);
        aabbTree = tree;
        aabbTree->accelerate_distance_queries();
}

Surface_mesh& Mymesh::get_raw_mesh()
{
    return mesh;
}

Tree* Mymesh::get_aabb_tree()
{
    return aabbTree.get();
}

Mymesh::Mymesh(const std::string &file_path)
{
    load_from_off(file_path);
    triangulate_mesh();
}

Mymesh::Mymesh(const Surface_mesh sm)
{
    this->mesh = sm;
    triangulate_mesh();    
}

Mymesh::Mymesh(double center_x, double center_y, double center_z, double dimension_x, double dimension_y, double dimension_z)
{
    double min_x = center_x - dimension_x/2, min_y = center_y - dimension_y/2, min_z = center_z - dimension_z/2;
    double max_x = center_x + dimension_x/2, max_y = center_y + dimension_y/2, max_z = center_z + dimension_z/2;

    Point v000(min_x, min_y, min_z);
    Point v100(max_x, min_y, min_z);
    Point v010(min_x, max_y, min_z);
    Point v001(min_x, min_y, max_z);
    Point v110(max_x, max_y, min_z);
    Point v101(max_x, min_y, max_z);
    Point v011(min_x, max_y, max_z);
    Point v111(max_x, max_y, max_z);

    std::vector<Point> vertices = {v000, v100, v110, v010, v001, v101, v111, v011};
    std::vector<vertex_descriptor> vd;

    Surface_mesh tissue_mesh;
    for (auto &p: vertices)
    {
        vertex_descriptor u = tissue_mesh.add_vertex(p);
        vd.push_back(u);
    } 

    tissue_mesh.add_face(vd[3], vd[2], vd[1], vd[0]);
    tissue_mesh.add_face(vd[4], vd[5], vd[6], vd[7]);
    tissue_mesh.add_face(vd[4], vd[7], vd[3], vd[0]);
    tissue_mesh.add_face(vd[1], vd[2], vd[6], vd[5]);
    tissue_mesh.add_face(vd[0], vd[1], vd[5], vd[4]);
    tissue_mesh.add_face(vd[2], vd[3], vd[7], vd[6]);

    this->mesh = tissue_mesh;
    triangulate_mesh();
}

std::string Mymesh::to_wkt() 
{
    std::stringstream ss;
	ss<<"POLYHEDRALSURFACE Z (";
    bool lfirst = true;

    for (face_descriptor f: mesh.faces())
    {
        if (lfirst) lfirst = false;
        else ss << ",";

        ss << "((";
        bool first = true;
        Point firstpoint;

        for (vertex_descriptor v: vertices_around_face(mesh.halfedge(f), mesh))
        {
            if (first) 
            {
                firstpoint = mesh.point(v);
                first = false;
            }
            else ss << ",";
            Point p = mesh.point(v);
            ss << p[0] << " " << p[1] << " " << p[2];
        }

        ss << "," << firstpoint[0] << " " << firstpoint[1] << " " << firstpoint[2];
        ss << "))";

    }
    ss << ")";

    return ss.str();
}