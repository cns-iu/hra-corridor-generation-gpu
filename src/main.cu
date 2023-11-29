#include "corridor.cuh"
#include <stdlib.h>
#include <stdio.h>

// global variables
std::unordered_map<std::string, Organ> total_body;

int main(int argc, char **argv)
{
    // print device information
    int count;
    cudaDeviceProp prop;
    cudaGetDeviceCount(&count);

    for (int i = 0; i < count; i++)
    {
        cudaGetDeviceProperties( &prop, i);
        printf( " --- General Information for device %d ---\n", i ); 
        printf( "Name: %s\n", prop.name );

        printf( "Threads in warp: %d\n", prop.warpSize );
        printf( "Max threads per block: %d\n", prop.maxThreadsPerBlock );
        printf( "Max thread dimensions: (%d, %d, %d)\n", prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2] );
        printf( "Max grid dimensions: (%d, %d, %d)\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2] ); 
        
        printf( "\n" );
    }

    // load arguments
    std::string body_path = std::string(argv[1]);

    // load all organs
    loadAllOrganModels(body_path, total_body);

    // rui location: http://purl.org/ccf/1.5/64c5c92d-dc07-40c4-a702-cbef17753fa6
    std::vector<std::pair<int, float>> collision_detection_result;
    collision_detection_result.push_back(std::make_pair(3, 0.031));
    collision_detection_result.push_back(std::make_pair(13, 0.961));

    AATissue example_tissue(0, 0, 0, 0.002, 0.003, 0.002);
    test_corridor_for_multiple_AS(example_tissue, collision_detection_result, total_body["VH_F_Kidney_R"], 0.05);


    return 0;
}