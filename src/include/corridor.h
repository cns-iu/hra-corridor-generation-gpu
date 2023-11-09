#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>

#include <stdlib.h>
#include <stdio.h>


// Axis-Aligned Tissue
class AATissue{
    public:
        AATissue() = default;
        AATissue(float center_x, float center_y, float center_z, float dimension_x, float dimension_y, float dimension_z);

    public:
        float center_x, center_y, center_z;
        float dimension_x, dimension_y, dimension_z;

};

class Organ{
    public:
        // meshes_vector is the concatenation of a 3d points. Three 3d points form a face. 
        // (offset[i+1] - offset[i]) * 3 form the i-th mesh
        std::vector<float3> meshes_vector;
        std::vector<uint> offset;
};