# ifndef KIRCHHOFF_CUH
# define KIRCHHOFF_CUH

# include "migration.hpp"

# define PI 3.14159265359

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

__global__ void cross_correlation(float * Ts, float * Tr, float * image, float * seismic, float aperture_x, float aperture_y, float cmp_x, float cmp_y, int spread, int nx, int ny, int nz, int nt, float dt, float dx, float dy, float dz);

# endif