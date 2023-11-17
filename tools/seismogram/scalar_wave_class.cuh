# ifndef ACOUSTIC_CUH
# define ACOUSTIC_CUH

# include <chrono>
# include <cuda_runtime.h>
# include <sys/resource.h>

# include "../../src/io/io.hpp"

class Acoustic
{
    std::chrono::_V2::system_clock::time_point ti, tf;

    int wavefieldBlocks, seismogramBlocks;
    int output_samples, threadsPerBlock;
    int time_id, source_id, total_nodes;
    int nt, nx, ny, nz, nb, model_size;
    int nxx, nyy, nzz, volume_size; 

    int RAM, vRAM, ivRAM;
    float dx, dy, dz, dt;

    // CPU array only

    int * sx = nullptr;
    int * sy = nullptr;
    int * sz = nullptr;
    
    float * vp = nullptr;
    float * dtvp2 = nullptr;

    float * output_data = nullptr;

    // GPU array only

    int * gx = nullptr;
    int * gy = nullptr;
    int * gz = nullptr;

    float * U = nullptr;
    float * Unew = nullptr;
    float * Uold = nullptr;
    float * dtVp2 = nullptr;
    
    float * damp1D = nullptr;
    float * damp2D = nullptr;
    float * damp3D = nullptr;

    float * wavelet = nullptr;
    float * seismogram = nullptr;

    void get_RAM_usage();
    void get_GPU_usage();
    void get_GPU_initMem();

    void import_parameters();

    void prepare_wavelet();
    void prepare_dampers();
    void prepare_velocity();
    void prepare_geometry();    
    void prepare_wavefield();    

    void expand_model_boundary();

public:

    int shot_id;

    int total_shots;

    std::string file;

    void set_patameters();

    void info_message();
    void initial_setup();
    void forward_solver();

    void set_runtime();
    void get_runtime();

    void export_outputs();

    void free_space();   
}; 

__global__ void update_wavefield(float * Unew, float * U, float * Uold, int volume_size);
__global__ void compute_wavefield(float * Unew, float * U, float * Uold, float * dtVp2, float * damp1D, float * damp2D, float * damp3D, float * wavelet, int source_id, int time_id, float dx, float dy, float dz, int nxx, int nyy, int nzz, int nb);
__global__ void compute_seismogram(float * Unew, float * seismogram, int * gx, int * gy, int * gz, int total_nodes, int time_id, int nt, int nxx, int nzz);
__device__ float get_boundary_damper(float * damp1D, float * damp2D, float * damp3D, int i, int j, int k, int nxx, int nyy, int nzz, int nb);

# endif