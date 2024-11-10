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

__global__ void cross_correlation(float * seismic, float * Ts, float * Tr, float * image, int nPoints, int spread, int nt, float dt);

# endif