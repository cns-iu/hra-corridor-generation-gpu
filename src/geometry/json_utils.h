#pragma once

#include <iostream>
#include <map>
#include <unordered_map>
#include <set>
#include <string>

#include<boost/filesystem.hpp>

#include <cpprest/http_listener.h>
#include <cpprest/json.h>
// #pragma comment(lib, "cpprest_2_10")


using namespace web;
using namespace web::http;
// using namespace web::http::experimental::listener;

struct Placement {
    Placement() = default;
    Placement(const std::string &src, const std::string &tgt, double x_scl, double y_scl, double z_scl, 
            double x_trans, double y_trans, double z_trans, 
            double x_rot, double y_rot, double z_rot): source(src), target(tgt), x_scaling(x_scl), y_scaling(y_scl), z_scaling(z_scl), x_translation(x_trans), y_translation(y_trans), z_translation(z_trans) {};

    std::string source;
    std::string target;
    double x_scaling, y_scaling, z_scaling;
    double x_translation, y_translation, z_translation;
    double x_rotation, y_rotation, z_rotation;

};


// placement transformation from old version models to new version models
void load_organ_transformation(const std::string &json_file_path, std::unordered_map<std::string, Placement> &mapping_placement);

json::value read_json_file(std::string const &json_file_path);

// placement structure, V1.2 -> V1.1 -> V1.0
void adjust_placement_target_ancestor(std::unordered_map<std::string, Placement> &mapping_placement);

void print_mapping_placement(std::unordered_map<std::string, Placement> &mapping_placement);