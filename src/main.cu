// #include "corridor.cuh"
#include <stdlib.h>
#include <stdio.h>


int main()
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

    return 0;
}