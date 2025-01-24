# ifndef KIRCHHOFF_CUH
# define KIRCHHOFF_CUH

# include "migration.hpp"

class Kirchhoff : public Migration
{
private:

    int nBlocks;
    int nThreads;

    float * d_Tr = nullptr;
    float * d_Ts = nullptr;

    float * d_image = nullptr;
    
    float * d_seismic = nullptr;

    void set_specifications();
    void run_cross_correlation();
};

__global__ void cross_correlation(float * Ts, float * Tr, float * image, float * seismic, float aperture_x, float aperture_y, 
                                  float cmp_x, float cmp_y, int spread, int nx, int ny, int nz, float dx, float dy, float dz, 
                                  int snx, int sny, int snz, float sdx, float sdy, float sdz, float scale, int nt, float dt);

__device__ float cubic1d(float P[4], float dx);
__device__ float cubic2d(float P[4][4], float dx, float dy);
__device__ float cubic3d(float P[4][4][4], float dx, float dy, float dz);

# endif