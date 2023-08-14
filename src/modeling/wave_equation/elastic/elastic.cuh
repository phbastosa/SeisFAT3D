# ifndef ELASTIC_CUH
# define ELASTIC_CUH

# include "../wave.cuh"

class Elastic : public Wave
{
protected:

    float * B = nullptr;
    float * L = nullptr;
    float * M = nullptr;

    float * Vx = nullptr;
    float * Vy = nullptr;
    float * Vz = nullptr;

    float * Txx = nullptr;
    float * Tyy = nullptr;
    float * Tzz = nullptr;
    float * Txz = nullptr;
    float * Tyz = nullptr;
    float * Txy = nullptr;
    
    void set_wavelet();
    void set_shear_model();
    void set_density_model();
    void set_model_boundaries();
    void set_modeling_volumes();
public:

    float * Vs;
    float * Rho;

    void set_parameters();
    void info_message();
    void forward_solver();
    void initial_setup();
    void free_space();
};

__global__ void compute_velocity(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * B, float * wavelet, int source_id, int time_id, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz);
__global__ void compute_stress(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * Pressure, float * M, float * L, float * damp1D, float * damp2D, float * damp3D, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz, int nb, int nbzu);
__device__ float get_boundary_damper(float * damp1D, float * damp2D, float * damp3D, int i, int j, int k, int nxx, int nyy, int nzz, int nb, int nbzu);

# endif