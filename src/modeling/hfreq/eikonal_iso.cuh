# ifndef EIKONAL_ISO_CUH
# define EIKONAL_ISO_CUH

# include "eikonal.hpp"

class Eikonal_Iso : public Eikonal
{
private:

    float dx2i, dy2i, dz2i, dsum;
    float dz2dx2, dz2dy2, dx2dy2;

    int total_levels;
    int nSweeps, meshDim;
    int nThreads, nBlocks;

    float * d_T = nullptr;
    float * d_S = nullptr;

    int * d_sgnv = nullptr;
    int * d_sgnt = nullptr;

    void set_properties();
    void set_conditions();

    int iDivUp(int a, int b);

public:

    void forward_solver();
};

__global__ void inner_sweep(float * S, float * T, int * sgnt, int * sgnv, int sgni, int sgnj, int sgnk, 
                            int level, int xOffset, int yOffset, int xSweepOffset, int ySweepOffset, int zSweepOffset, 
                            int nxx, int nyy, int nzz, float dx, float dy, float dz, float dx2i, float dy2i, float dz2i, 
                            float dz2dx2, float dz2dy2, float dx2dy2, float dsum);

# endif