# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include "../modeling.hpp"

class Eikonal : public Modeling
{
private:

    void set_boundaries();
    void set_specifications();

    float cubic1d(float P[4], float dx);
    float cubic2d(float P[4][4], float dx, float dy);
    float cubic3d(float P[4][4][4], float dx, float dy, float dz);

protected:

    int total_levels;
    int nSweeps, meshDim;
    int nThreads, nBlocks;

    float dx2i, dy2i, dz2i, dsum;
    float dz2dx2, dz2dy2, dx2dy2;

    int * d_sgnv = nullptr;
    int * d_sgnt = nullptr;

    float * d_T = nullptr;
    float * d_S = nullptr;

    void compute_seismogram();

    virtual void set_conditions() = 0;
    virtual void set_properties() = 0;

    int iDivUp(int a, int b);

public:

    void propagation();
    void initialization();
    
    virtual void forward_solver() = 0;

    void export_synthetic_data();
};

__global__ void inner_sweep(float * S, float * T, int * sgnt, int * sgnv, int sgni, int sgnj, int sgnk, 
                            int level, int xOffset, int yOffset, int xSweepOffset, int ySweepOffset, int zSweepOffset, 
                            int nxx, int nyy, int nzz, float dx, float dy, float dz, float dx2i, float dy2i, float dz2i, 
                            float dz2dx2, float dz2dy2, float dx2dy2, float dsum);

# endif