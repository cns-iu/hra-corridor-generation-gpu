#include "utils.h"

namespace fs = boost::filesystem;
#define PI 3.14159265

void tissue_transform(std::unordered_map<std::string, double> &params, Surface_mesh &tissue_mesh, std::vector<Point> &points, int resolution=10)
{
    double x_dimension = params["x_dimension"];
    double y_dimension = params["y_dimension"];
    double z_dimension = params["z_dimension"];

    double x_translation = params["x_translation"];
    double y_translation = params["y_translation"];
    double z_translation = params["z_translation"];
    
    double x_scaling = params["x_scaling"];
    double y_scaling = params["y_scaling"];
    double z_scaling = params["z_scaling"];

    double x_rotation = params["x_rotation"];
    double y_rotation = params["y_rotation"];
    double z_rotation = params["z_rotation"];

    double x_origin = params["x_origin"];
    double y_origin = params["y_origin"];
    double z_origin = params["z_origin"];

    Eigen::Vector3d origin(x_origin, y_origin, z_origin);

    double min_x = -x_dimension/2.0 * x_scaling, min_y = -y_dimension/2.0 * y_scaling, min_z = -z_dimension/2.0 * z_scaling;
    double max_x = x_dimension/2.0 * x_scaling, max_y = y_dimension/2.0 * y_scaling, max_z = z_dimension/2.0 * z_scaling;

    Eigen::Vector3d v000(min_x, min_y, min_z);
    Eigen::Vector3d v100(max_x, min_y, min_z);
    Eigen::Vector3d v010(min_x, max_y, min_z);
    Eigen::Vector3d v001(min_x, min_y, max_z);
    Eigen::Vector3d v110(max_x, max_y, min_z);
    Eigen::Vector3d v101(max_x, min_y, max_z);
    Eigen::Vector3d v011(min_x, max_y, max_z);
    Eigen::Vector3d v111(max_x, max_y, max_z);
    
    std::vector<Eigen::Vector3d> vertices = {v000, v100, v110, v010, v001, v101, v111, v011};

    Eigen::MatrixXd R_x(3, 3);
    Eigen::MatrixXd R_y(3, 3);
    Eigen::MatrixXd R_z(3, 3);
    // Eigen::MatrixXd z_refect(3, 3);

    double radius_x = x_rotation * PI / 180;
    double radius_y = y_rotation * PI / 180;
    double radius_z = z_rotation * PI / 180;

    R_x << 1.0, 0.0, 0.0, 
        0.0, cos(radius_x), -sin(radius_x),
        0.0, sin(radius_x), cos(radius_x);
    R_y << cos(radius_y), 0.0, sin(radius_y),
        0.0, 1.0, 0.0,
        -sin(radius_y), 0.0, cos(radius_y);
    R_z << cos(radius_z), -sin(radius_z), 0.0, 
        sin(radius_z), cos(radius_z), 0.0,
        0.0, 0.0, 1.0;


    // z_refect << 1, 0, 0, 0, 1, 0, 0, 0, -1;
    
    Eigen::MatrixXd R = R_x * R_y * R_z;
    Eigen::Vector3d T(x_translation, y_translation, z_translation);

    // std::cout << "translation vector: " << T << " " << "origin" << origin << std::endl;
   
    std::vector<vertex_descriptor> vd;
    
    for (int i = 0; i < vertices.size(); i++)
    {   
        vertices[i] = (R * vertices[i] + T)/1000.0 + origin;
    }

    for (auto &vec: vertices)
    {
        Point p(vec(0), vec(1), vec(2));
        vertex_descriptor u = tissue_mesh.add_vertex(p);
        vd.push_back(u);
    }

    tissue_mesh.add_face(vd[3], vd[2], vd[1], vd[0]);
    tissue_mesh.add_face(vd[4], vd[5], vd[6], vd[7]);
    tissue_mesh.add_face(vd[4], vd[7], vd[3], vd[0]);
    tissue_mesh.add_face(vd[1], vd[2], vd[6], vd[5]);
    tissue_mesh.add_face(vd[0], vd[1], vd[5], vd[4]);
    tissue_mesh.add_face(vd[2], vd[3], vd[7], vd[6]);

    std::ofstream tissue_mesh_off("tissue_mesh.off");
    tissue_mesh_off << tissue_mesh;
    tissue_mesh_off.close();

    double delta_x = (max_x - min_x) / resolution, delta_y = (max_y - min_y) / resolution, delta_z = (max_z - min_z) / resolution;    
    double center_x, center_y, center_z;

    for (int i = 0; i < resolution; i++)
        for (int j = 0; j < resolution; j++)
            for (int k = 0; k < resolution; k++)
            {
                center_x = min_x + (i + 0.5) * delta_x;
                center_y = min_y + (j + 0.5) * delta_y;
                center_z = min_z + (k + 0.5) * delta_z; 
                Eigen::Vector3d vec(center_x, center_y, center_z);
                // Affine transformation
                vec = (R*vec + T)/1000.0 + origin;
                Point p(vec(0), vec(1), vec(2));
                points.push_back(p);
            }

}

void load_all_organs(const std::string &body_path, std::unordered_map<std::string, std::vector<Mymesh>> &total_body)
{
    // std::string body_path = "/home/luchen/collision_detection_kidney_2/models/AS";
    std::ofstream output_wkt;
    output_wkt.open("output.wkt");

    for (fs::directory_entry& organ_path : fs::directory_iterator(body_path)) 
    {
        std::string organ_name = organ_path.path().stem().string();
        std::cout << organ_name << std::endl;   
        for (fs::directory_entry& AS : fs::directory_iterator(organ_path)) 
        {
            std::string file_path = AS.path().string();
            total_body[organ_name].push_back(Mymesh(file_path));
        }

        for (auto &AS: total_body[organ_name]) {
            AS.create_aabb_tree();
            // std::string as_wkt = AS.to_wkt();
            // output_wkt << organ_name << "|" << AS.label << "|" << as_wkt << "\n";

        }
    }

    output_wkt.close();

}


void gen_origin(const std::string &organ_origins_file, std::unordered_map<std::string, Eigen::Vector3d> &organ_origins)
{

    Eigen::Vector3d VHFRightKidney_origin(-0.1087215, 0.1985958, -0.1344381);
    Eigen::Vector3d VHFLeftKidney_origin(0.04005157, 0.1757761, -0.1393864);
    Eigen::Vector3d VHMRightKidney_origin(-0.09692713, 0.22901, -0.05137664);
    Eigen::Vector3d VHMLeftKidney_origin(0.03373008, 0.2498251, -0.06865281);

    organ_origins["#VHFRightKidney"] = VHFRightKidney_origin;
    organ_origins["#VHFLeftKidney"] = VHFLeftKidney_origin;
    organ_origins["#VHMRightKidney"] = VHMRightKidney_origin;
    organ_origins["#VHMLeftKidney"] = VHMLeftKidney_origin;

    std::ifstream origins(organ_origins_file);

    if (!origins.is_open()) throw std::runtime_error("could not open " + organ_origins_file);
    
    std::string line;
    if (origins.good())
    {
        std::getline(origins, line); //skip the first line

        std::vector<std::string> row;
        std::string word;

        while (std::getline(origins, line))
        {
            row.clear();
            std::stringstream ss(line);
            while (std::getline(ss, word, ','))
            {
                row.push_back(word);
            }

            auto target = row[0];
            auto x_translation = std::stod(row[1]);
            auto y_translation = std::stod(row[2]);
            auto z_translation = std::stod(row[3]);

            // handle suffix
            // size_t len = target.size();
            // if (target[len - 4] == 'V' && target[len-2] == '.') target = target.substr(0, len-4);

            Eigen::Vector3d origin(x_translation, y_translation, z_translation);
            organ_origins[target] = origin;


        }

    }

}

void load_ASCT_B(const std::string &file_path, std::unordered_map<std::string, std::string> &mapping, std::unordered_map<std::string, SpatialEntity> &mapping_node_spatial_entity)
{
    
    mapping["#VHFRightKidney"] = "VH_F_Kidney_R";
    mapping["#VHFLeftKidney"] = "VH_F_Kidney_L";
    mapping["#VHMRightKidney"] = "VH_M_Kidney_R";
    mapping["#VHMLeftKidney"] = "VH_M_Kidney_L";
    

    std::ifstream asct_b(file_path);

    if (!asct_b.is_open()) throw std::runtime_error("could not open asct_b table!");

    std::string line, colname;
    if (asct_b.good())
    {
        for (int i = 0; i < 10; i++) std::getline(asct_b, line);

        std::getline(asct_b, line);
        std::stringstream ss(line);

        while (std::getline(ss, colname, ','));
    }

    std::vector<std::string> row;
    std::string word;
    while (std::getline(asct_b, line))
    {
        row.clear();
        std::stringstream ss(line);
        while (std::getline(ss, word, ','))
        {
            row.push_back(word);
        }

        auto anatomical_structure_of = row[0];
        auto source_spatial_entity = row[1];
        auto node_name = row[2];
        auto label = row[3];
        auto ontologyID = row[4];
        auto representation_of = row[5];
        auto glb_file = row[6];

        if (anatomical_structure_of != "-")
        {
            size_t n = glb_file.size();
            char last_char = glb_file[n - 1];
            if (last_char == '\r' || last_char == '\n' || last_char == '\t') glb_file = glb_file.substr(0, n - 1);

            // end with V1.1 or V1.2 or any V*.* 
            // if (n > 4 && anatomical_structure_of[n-4] == 'V' && anatomical_structure_of[n-2] == '.') 
            //     anatomical_structure_of = anatomical_structure_of.substr(0, n-4);
            
            mapping[anatomical_structure_of] = glb_file;
        }
        
        
        if (node_name != "-")
        {
            // mapping_node_ontology[node_name] = representation_of;
            // mapping_node_label[node_name] = label;
            SpatialEntity spatial_entity(anatomical_structure_of, source_spatial_entity, node_name, label, representation_of, glb_file);
            mapping_node_spatial_entity[node_name] = spatial_entity; 

        }


    }


    

}

std::string organ_split(const std::string &url)
{

    size_t len = url.size();
    
    // end with V1.1 or V1.2 or any V*.* 
    // if (url[len-4] == 'V' && url[len-2] == '.') tmp = url.substr(0, len-4);
    // else tmp = url;

    // start from # (including #)
    size_t start = url.find("#"); 
    return url.substr(start);

}

void output_corridor(Surface_mesh &mesh, std::string rui_location_id, std::string output_corridor_dir)
{
    if (!fs::exists(output_corridor_dir)) fs::create_directory(output_corridor_dir);
    size_t start = rui_location_id.find_last_of('/') + 1;
    std::string abs_file_path = output_corridor_dir + '/' + rui_location_id.substr(start) + ".off";
    CGAL::IO::write_polygon_mesh(abs_file_path, mesh, CGAL::parameters::stream_precision(17));            

}


std::string output_corridor_glb(Surface_mesh &mesh, std::string rui_location_id, std::string output_corridor_dir)
{
    if (!fs::exists(output_corridor_dir)) fs::create_directory(output_corridor_dir);
    size_t start = rui_location_id.find_last_of('/') + 1;
    std::string off_file_path = output_corridor_dir + '/' + rui_location_id.substr(start) + ".off";
    CGAL::IO::write_polygon_mesh(off_file_path, mesh, CGAL::parameters::stream_precision(17));

    std::string glb_file_path = output_corridor_dir + '/' + rui_location_id.substr(start) + ".glb";
    // execute python script for the format conversion
    std::string command = "python " + off_file_path + " " + glb_file_path + " ./scripts/glb_converter.py";
    system(command.c_str());

    // generate binary string
    std::stringstream buffer;
    if (fs::exists(glb_file_path))
    {
        std::ifstream if_glb(glb_file_path);
        buffer << if_glb.rdbuf();
    }       

    return buffer.str();

}

std::string output_corridor_glb(Surface_mesh &corridor_mesh, std::string rui_location_id)
{

    std::string output_corridor_dir = "./tmp_corridors";
    if (!fs::exists(output_corridor_dir)) fs::create_directory(output_corridor_dir);
    size_t start = rui_location_id.find_last_of('/') + 1;
    std::string off_file_path = output_corridor_dir + '/' + rui_location_id.substr(start) + ".off";
    std::ofstream corridor_output(off_file_path);
    corridor_output << corridor_mesh;

     // Create an Assimp importer
    Assimp::Importer importer;
    // Read the OFF file
    const aiScene* scene = importer.ReadFile(off_file_path, aiProcess_Triangulate); //| aiProcess_FlipUVs);
    std::stringstream buffer;
    
    if (scene) {
        // Create an Assimp exporter
        Assimp::Exporter exporter;
        std::string glb_file_path = output_corridor_dir + '/' + rui_location_id.substr(start) + ".glb";
        // Export the scene to GLB format
        exporter.Export(scene, "glb2", glb_file_path); //aiProcess_Triangulate | aiProcess_FlipUVs);

        // generate binary string
        if (fs::exists(glb_file_path))
        {
            std::ifstream if_glb(glb_file_path);
            buffer << if_glb.rdbuf();
        }

        std::cout << "Conversion successful. GLB file saved to: " << glb_file_path << std::endl;
    } else {
        std::cerr << "Error loading the OFF file: " << importer.GetErrorString() << std::endl;
    }       

    return buffer.str();

}

void comparison_CPU_GPU(std::vector<Point> &CPU_points, std::vector<Point> &GPU_points)
{
    if (CPU_points.size() != GPU_points.size()) 
    {
        std::cout << "[size] GPU computation is not the same with CPU! CPU size: " << CPU_points.size() << " GPU size: " << GPU_points.size() << std::endl;
    }
    // for (int i = 0; i < CPU_points.size(); i++)
    // {
    //     auto &p1 = CPU_points[i];
    //     auto &p2 = GPU_points[i];

    //     if (std::abs(p1[0] - p2[0]) > 1e-6 || std::abs(p1[1] - p2[1]) > 1e-6 || std::abs(p1[2] - p2[2]) > 1e-6)
    //     {
    //         std::cout << "[points] GPU computation is not the same with CPU!" << std::endl;
    //     } 
    // }
}
