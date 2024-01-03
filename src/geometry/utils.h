#pragma once
#include <sstream>
#include <iostream>
#include<fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <Eigen/Dense>
#include <math.h>
#include<boost/filesystem.hpp>

#include "mymesh.h"


struct SpatialEntity {

    SpatialEntity() = default;
    SpatialEntity(const std::string &aso, const std::string &sse, const std::string &nn, const std::string &lb, const std::string &ro, const std::string &gf): anatomical_structure_of(aso), source_spatial_entity(sse), node_name(nn), label(lb), representation_of(ro), glb_file(gf) {};

    std::string anatomical_structure_of;
    std::string source_spatial_entity;
    std::string node_name;
    std::string label;
    std::string representation_of;
    std::string glb_file;


};

//origins, meshes
void load_all_organs(const std::string &body_path, std::unordered_map<std::string, std::vector<Mymesh>> &total_body);

// generate origins from the whole model without any hole-filling.
void gen_origin(const std::string &whole_model_path, std::unordered_map<std::string, Eigen::Vector3d> &organ_origins);

//including x_scaling, x_rotation, x_translation, x_origin 
void tissue_transform(std::unordered_map<std::string, double> &params, Surface_mesh &tissue_mesh, std::vector<Point> &points, int resolution);

std::string organ_split(const std::string &url);

void load_ASCT_B(const std::string &file_path, std::unordered_map<std::string, std::string> &mapping, std::unordered_map<std::string, SpatialEntity> &mapping_node_spatial_entity);

void output_corridor(Surface_mesh &mesh, std::string rui_location_id, std::string output_corridor_dir);