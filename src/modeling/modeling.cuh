# ifndef MODELING_HPP
# define MODELING_HPP

# include <chrono>
# include <cuda_runtime.h>
# include <sys/resource.h>

# include "../geometry/geometry.hpp"
# include "../geometry/regular/regular.hpp"
# include "../geometry/circular/circular.hpp"

class Modeling
{
private:

    int RAM, vRAM, ivRAM;
    void get_RAM_usage();
    void get_GPU_usage();
    void get_GPU_initMem();

    std::chrono::_V2::system_clock::time_point ti, tf;

    void check_geometry_overflow();

protected:

    float t0;    
    
    float dx, dy, dz;
    int nx, ny, nz, nt; 
    int nPoints, volsize;
    int nit, nxx, nyy, nzz;

    float * S = nullptr;
    float * V = nullptr;
    float * T = nullptr;

    float * d_T = nullptr;
    float * d_S = nullptr;

    int * d_sgnv = nullptr;
    int * d_sgnt = nullptr;

    int threadsPerBlock;
    int nSweeps, meshDim;
    int totalLevels, padb;
    
    float dx2i, dy2i, dz2i, dsum; 
    float dx2dy2, dz2dx2, dz2dy2;

    int receiver_output_samples;
    int wavefield_output_samples;

    bool export_receiver_output;
    bool export_wavefield_output;

    float * receiver_output = nullptr;
    float * wavefield_output = nullptr;

    std::string receiver_output_file;
    std::string wavefield_output_file;
    std::string receiver_output_folder;
    std::string wavefield_output_folder;

    Geometry * geometry;

    void get_travelTimes();
    void get_firstArrivals();

    int iDivUp(int a, int b); 

public: 

    int shot_id;
    int total_shots;
    int total_nodes;
    
    std::string file;

    void set_runtime();
    void get_runtime();
    void info_message();
    void export_outputs();

    void set_parameters(); 
    void forward_solver();
    void initial_setup();
    void build_outputs();
    void set_slowness();
    void free_space();
};

__global__ void fast_sweeping_kernel(float * S, float * T, int * sgnt, int * sgnv, int sgni, int sgnj, int sgnk, 
                                     int level, int xOffset, int yOffset, int xSweepOffset, int ySweepOffset, int zSweepOffset, 
                                     int nxx, int nyy, int nzz, float dx, float dy, float dz, float dx2i, float dy2i, float dz2i, 
                                     float dz2dx2, float dz2dy2, float dx2dy2, float dsum);

# endif