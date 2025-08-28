# ifndef MODELING_CUH
# define MODELING_CUH

# include <cuda_runtime.h>

# include "../geometry/geometry.hpp"

# define NSWEEPS 8
# define MESHDIM 3

# define COMPRESS 65535

typedef unsigned short int uintc; 

class Modeling
{
private:
    
    void set_eikonal();
    void set_properties();

    int iDivUp(int a, int b);

    float cubic1d(float P[4], float dx);
    float cubic2d(float P[4][4], float dx, float dy);
    float cubic3d(float P[4][4][4], float dx, float dy, float dz);

protected:

    int total_levels;
    int nThreads, nBlocks;

    float dx2i, dy2i, dz2i, dsum;
    float dz2dx2, dz2dy2, dx2dy2;

    int * d_sgnv = nullptr;
    int * d_sgnt = nullptr;

    virtual void set_conditions() = 0;
    
    void compression(float * input, uintc * output, int volsize, float &max_value, float &min_value);

public:

    float dx, dy, dz;
    int nxx, nyy, nzz, volsize;
    int nx, ny, nz, nb, nPoints;
    int srcId, recId, sIdx, sIdy, sIdz;

    float sx, sy, sz;

    float * S = nullptr;
    float * T = nullptr;

    float * d_T = nullptr;
    float * d_S = nullptr;

    float * seismogram = nullptr;

    int max_spread;
    Geometry * geometry;
    
    std::string parameters;
    std::string data_folder;
    std::string modeling_type;
    std::string modeling_name;

    void set_parameters();
    void initialization();
    void eikonal_solver();
    void set_shot_point();
    void show_information();    
    void compute_seismogram();

    void copy_slowness_to_device();

    void expand_boundary(float * input, float * output);
    void reduce_boundary(float * input, float * output);
    
    virtual void time_propagation() = 0;

    void export_seismogram();
};

__global__ void time_set(float * T, int volsize);

__global__ void time_init(float * T, float * S, float sx, float sy, float sz, float dx, float dy, 
                          float dz, int sIdx, int sIdy, int sIdz, int nxx, int nzz, int nb);

__global__ void inner_sweep(float * S, float * T, int * sgnt, int * sgnv, int sgni, int sgnj, int sgnk, 
                            int level, int xOffset, int yOffset, int xSweepOffset, int ySweepOffset, int zSweepOffset, 
                            int nxx, int nyy, int nzz, float dx, float dy, float dz, float dx2i, float dy2i, float dz2i, 
                            float dz2dx2, float dz2dy2, float dx2dy2, float dsum);

# endif