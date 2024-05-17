#include "algo.h"
#include "utils.h"
#include "json_utils.h"
#include "corridor_cpu.h"
#include "corridor.cuh"
#include "aabbtree_gpu.h"
#include "corridor_with_spatial_index.cuh"
#include <chrono>

#include <cpprest/http_listener.h>
#include <cpprest/json.h>
#pragma comment(lib, "cpprest_2_10")

using namespace web;
using namespace web::http;
using namespace web::http::experimental::listener;

#include <iostream>
#include <map>
#include <unordered_map>
#include <set>
#include <string>

//global variables
std::map<utility::string_t, utility::string_t> dictionary;
std::unordered_map<std::string, Eigen::Vector3d> organ_origins;                     //origins of organs
std::unordered_map<std::string, std::string> mapping;                               //mapping from standard organ name(e.g., #VHFLeftKidney) to glb file name without suffix(e.g., VH_F_Kidney_L)
std::unordered_map<std::string, std::vector<Mymesh>> total_body;                    //mapping from organ name(glb file name) to vector of meshes of a certain organ
std::unordered_map<std::string, SpatialEntity> mapping_node_spatial_entity;         // mapping from AS to its information in asct-b table 
std::unordered_map<std::string, Placement> mapping_placement;                       // mapping from source plcement to target placement(including rotation, scaling, translation parameters)
std::unordered_map<std::string, Organ> total_body_gpu;                              // gpu version of total_body: key is organ name, value is Organ object
std::unordered_map<std::string, std::vector<myspatial::AABBTree*>> total_body_gpu_with_spatial_index; // gpu version of total_body with spatial index. key is the organ name, value is the pointer of aabbtrees

// Point format convertion from gpu (float3) to cpu (Point in CGAL)
std::vector<Point> convert_point_gpu_to_cpu(ResultContainer &result_container, AATissue &example_tissue)
{
   std::vector<Point> point_cloud;
   float3* p_corridor_arr = result_container.corridor_array;
   bool* p_point_is_in_corridor_array = result_container.point_is_in_corridor_array;
   
   float example_d_x = example_tissue.dimension_x;
   float example_d_y = example_tissue.dimension_y;
   float example_d_z = example_tissue.dimension_z;

   float c_x, c_y, c_z;
   for (int i = 0; i < 40*40*40; i++) {
      if (p_point_is_in_corridor_array[i])
      {
         float3 p = p_corridor_arr[i];
         c_x = p.x;
         c_y = p.y;
         c_z = p.z;
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


//display json
void display_json(
   json::value const & jvalue,
   utility::string_t const & prefix)
{
   // std::cout << prefix << jvalue.serialize() << std::endl;
}

//parse json
void parse_json(json::value const &jvalue, json::value &answer)
{
   // for test
   // std::string tissue_output_path = "tissue_mesh_2.off";
   try 
   {
         auto placement = jvalue.at("placement");
         std::unordered_map<std::string, double> params;
         auto target = placement.at("target").as_string();

         // for test
         // if (jvalue.has_field("tissue_block_id")) tissue_output_path = jvalue.at("tissue_block_id").as_string() + ".off";

         //extract parameters from json request
         params["x_dimension"] = jvalue.at("x_dimension").as_double();
         params["y_dimension"] = jvalue.at("y_dimension").as_double();
         params["z_dimension"] = jvalue.at("z_dimension").as_double();
         params["x_scaling"] = placement.at("x_scaling").as_double();
         params["y_scaling"] = placement.at("y_scaling").as_double();
         params["z_scaling"] = placement.at("z_scaling").as_double();
         params["x_translation"] = placement.at("x_translation").as_double();
         params["y_translation"] = placement.at("y_translation").as_double();
         params["z_translation"] = placement.at("z_translation").as_double();
         params["x_rotation"] = placement.at("x_rotation").as_double();
         params["y_rotation"] = placement.at("y_rotation").as_double();
         params["z_rotation"] = placement.at("z_rotation").as_double();

         if (!(mapping_placement.find(target) == mapping_placement.end()))
         {
            Placement &placement = mapping_placement[target];
            target = placement.target;
            params["x_translation"] *= placement.x_scaling;
            params["y_translation"] *= placement.y_scaling;
            params["z_translation"] *= placement.z_scaling;
            // other transformations here. 
         }

         auto reference_organ_name = organ_split(target); 

         // test for all organs
         if (mapping.find(reference_organ_name) == mapping.end()) 
         {
            std::cout << reference_organ_name << " doesn't exist in ASCT-B table!" << std::endl;
            return;
         }
         
         std::string organ_file_name = mapping[organ_split(target)];
         std::cout << "target url: " << target << " target: " << reference_organ_name << " " << "organ file name: " << organ_file_name << std::endl;

         Eigen::Vector3d origin = organ_origins[reference_organ_name];
         params["x_origin"] = origin(0);
         params["y_origin"] = origin(1);
         params["z_origin"] = origin(2);

         
         Surface_mesh tissue_mesh;
         std::vector<Point> points; //center of voxels inside the tissue block
         tissue_transform(params, tissue_mesh, points, 10);

         Mymesh my_tissue(tissue_mesh);

         // for test
         // std::ofstream tissue_mesh_off("./tissue_blocks/" + tissue_output_path);
         // tissue_mesh_off << tissue_mesh;
         // tissue_mesh_off.close();


         my_tissue.create_aabb_tree();

         //core function
         auto t1 = std::chrono::high_resolution_clock::now();
         std::vector<std::pair<int, double>> result = collision_detection_single_tissue_2(total_body[organ_file_name], my_tissue, points);
         auto t2 = std::chrono::high_resolution_clock::now();

         std::chrono::duration<double> duration2 = t2 - t1;
         // std::cout << "collision detection function running time: " << duration2.count() << " seconds" << std::endl;  

         auto result_bb = collision_detection_bb(total_body[organ_file_name], my_tissue);

         auto &target_organ = total_body[organ_file_name];
         //print result
         std::cout << "mesh collision detection result:\nlabel         percentage" << std::endl;
         for (auto s: result) {std::cout << target_organ[s.first].label << " " << s.second << std::endl;}
         std::cout << "bounding box collision detection result: \nlabel" << std::endl;
         for (auto s: result_bb) {std::cout << s << std::endl; }

         //construct the response
         double tissue_volume = params["x_dimension"] * params["y_dimension"] * params["z_dimension"];

         // GPU corridor generation based on collision detection result
         auto rui_location_id = jvalue.at("@id").as_string();
         std::cout << "rui location id: " << rui_location_id << std::endl;
         std::string output_corridor_dir = "./corridor_models";
         
         double tolerance = 0.05;

         Surface_mesh corridor_mesh;
         if (result.size() > 3) {
            corridor_mesh = my_tissue.get_raw_mesh();
            // glb = output_corridor_glb(my_tissue.get_raw_mesh(), rui_location_id, output_corridor_dir);
         }
         else if (result.size() == 2 || result.size() == 3)
         {           
            AATissue example_tissue_gpu(0, 0, 0, params["x_dimension"]/1000.0, params["y_dimension"]/1000.0, params["z_dimension"]/1000.0);
            Mytissue example_tissue_cpu(0.0, 0.0, 0.0, params["x_dimension"]/1000.0, params["y_dimension"]/1000.0, params["z_dimension"]/1000.0);
            // collision detection result: convert intersection volume (the second value), from double to float
            std::vector<std::pair<int, float>> float_collision_detection_result;
            for (auto s: result) float_collision_detection_result.push_back(std::make_pair(s.first, (float) s.second));
            
            // CPU baseline 1, using CGAL
            t1 = std::chrono::high_resolution_clock::now();
            auto baseline_points = create_point_cloud_corridor_for_multiple_AS(total_body[organ_file_name], example_tissue_cpu, result, tolerance);
            t2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_cpu = t2 - t1;
            std::cout << "CPU baseline: " << duration_cpu.count() << " seconds" << std::endl;

            // GPU naive implementation
            // note: loading total body (cpu and gpu) can merge together
            t1 = std::chrono::high_resolution_clock::now();
            ResultContainer result_container = test_corridor_for_multiple_AS(example_tissue_gpu, float_collision_detection_result, total_body_gpu[organ_file_name], tolerance);
            // from float3 arr to CGAL Point vector
            std::vector<Point> points = convert_point_gpu_to_cpu(result_container, example_tissue_gpu);
            t2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_gpu = t2 - t1;
            std::cout << "GPU corridor: " << duration_gpu.count() << " seconds" << std::endl;
            
            // GPU implementation with spatial index
            t1 = std::chrono::high_resolution_clock::now();
            ResultContainer result_container_2 = test_corridor_with_spatial_index_gpu(example_tissue_gpu, float_collision_detection_result, total_body_gpu_with_spatial_index[organ_file_name], tolerance);
            t2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration_gpu_spatial_index = t2 - t1;
            std::cout << "GPU corridor with spatial index: " << duration_gpu_spatial_index.count() << " seconds" << std::endl;

            // verification GPU results
            comparison_CPU_GPU(baseline_points, points);
            // using CGAL function to reconstruct mesh from points
            if (points.size() == 0) corridor_mesh = my_tissue.get_raw_mesh();
            else corridor_mesh = create_corridor_from_point_cloud(points);
            
            // construct the response of running time
            answer[U("CPU_time")] = json::value(duration_cpu.count());
            answer[U("GPU_time")] = json::value(duration_gpu.count());
            answer[U("GPU_time_spatial_index")] = json::value(duration_gpu_spatial_index.count());
         }
         else if (result.size() == 1)
         {
            corridor_mesh =  target_organ[result[0].first].get_raw_mesh();
         }
         
         std::string glb = output_corridor_glb(corridor_mesh, rui_location_id);

         // construct response for corridor
         answer[U("number_of_collisions")] = json::value(result.size());
         answer[U("corridor_glb")] = json::value::string(U(glb));
         answer[U("rui_location_id")] = json::value::string(U(rui_location_id));

   }
   catch(...)
   {
      std::cout << "catch exception in parse json function" << std::endl;
   }

}


void handle_get(http_request request)
{

   std::cout << "\nhandle GET" << std::endl;

   auto answer = json::value::object();

   for (auto const & p : dictionary)
   {
      answer[p.first] = json::value::string(p.second);
   }

   display_json(json::value::null(), "R: ");
   display_json(answer, "S: ");

   request.reply(status_codes::OK, answer);
}

void handle_request(http_request request, std::function<void(json::value const &, json::value &)> action)
{
   json::value answer;

   request
      .extract_json()
      .then([&answer, &action](pplx::task<json::value> task) {
         try
         {
            auto const & jvalue = task.get();
            display_json(jvalue, "Request: ");

            if (!jvalue.is_null())
            {
               action(jvalue, answer);
            }
            
            display_json(answer, "Response: ");
         }
         catch (http_exception const & e)
         {
            std::cout << e.what() << std::endl;
         }
      })
      .wait();

   // if (answer != json::value::null())
   //    request.reply(status_codes::OK, answer);
   // else
   //    request.reply(status_codes::OK, json::value::array());
   
   http_response response(status_codes::OK);
   response.headers().add(U("Access-Control-Allow-Origin"), U("*"));


   if (answer != json::value::null())
      response.set_body(answer);
   else
      response.set_body(json::value::array());
   
   request.reply(response);

}

void handle_post(http_request request)
{
   // TRACE("\nhandle POST\n");
   std::cout << "\nhandle POST" << std::endl;

   handle_request(
      request,
      parse_json
   );
}

void handle_options(http_request request)
{
   http_response response(status_codes::OK);
   response.headers().add(U("Allow"), U("GET, POST, OPTIONS"));
   response.headers().add(U("Access-Control-Allow-Origin"), U("*"));
   response.headers().add(U("Access-Control-Allow-Methods"), U("GET, POST, OPTIONS"));
   response.headers().add(U("Access-Control-Allow-Headers"), U("Content-Type"));
   request.reply(response);

}

int main(int argc, char **argv)
{

   // std::string organ_origins_file_path = "/home/catherine/data/model/organ_origins_meter.csv";
   // std::string asct_b_file_path = "/home/catherine/data/model/ASCT-B_3D_Models_Mapping.csv";
   // std::string body_path = "/home/catherine/data/model/plain_filling_hole";
   if (argc < 7)
   {
      std::cout << "Please provide the organ_origins_file_path, asct_b_file_path, body_path(model_path), reference_organ_json_file, server IP and port number!" << std::endl;
      return 0;
   }

   std::string organ_origins_file_path = std::string(argv[1]);
   std::string asct_b_file_path = std::string(argv[2]);
   std::string body_path = std::string(argv[3]);
   std::string reference_organ_json_file = std::string(argv[4]);
   std::string server_ip = std::string(argv[5]);
   std::string port = std::string(argv[6]);

   // load origins
   gen_origin(organ_origins_file_path, organ_origins);
   load_ASCT_B(asct_b_file_path, mapping, mapping_node_spatial_entity);
   load_all_organs(body_path, total_body);
   load_organ_transformation(reference_organ_json_file, mapping_placement);

   // load organ for GPU use, naive and spatial index
   loadAllOrganModels(body_path, total_body_gpu);
   loadAllOrganModels(body_path, total_body_gpu_with_spatial_index);

   http_listener listener("http://" + server_ip + ":" + port + "/get-collisions");

   listener.support(methods::GET,  handle_get);
   listener.support(methods::POST, handle_post);
   listener.support(methods::OPTIONS, handle_options);


   try
   {
      listener
         .open()
         .then([&listener]() {
            // TRACE("\nstarting to listen\n"); 
            std::cout << "\nstarting to listen" << std::endl;
            })
         .wait();

      while (true);
   }
   catch (std::exception const & e)
   {
      std::cout << e.what() << std::endl;
   }

   return 0;

}