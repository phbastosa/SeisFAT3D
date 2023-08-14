# ifndef ACOUSTIC_CUH
# define ACOUSTIC_CUH

# include "../wave.cuh" 

class Acoustic : public Wave
{
protected:

    float * K = nullptr;
    float * B = nullptr;
    float * Vx = nullptr; 
    float * Vy = nullptr; 
    float * Vz = nullptr; 

    void set_wavelet();
    void set_density_model();
    void set_model_boundaries();
    void set_modeling_volumes();

public:

    float * Rho = nullptr;

    void set_parameters();
    void info_message();
    void forward_solver();
    void initial_setup();
    void free_space();
};

__global__ void compute_velocity(float * Pressure, float * Vx, float * Vy, float * Vz, float * B, float * wavelet, int source_id, int time_id, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz);
__global__ void compute_pressure(float * Pressure, float * Vx, float * Vy, float * Vz, float * K, float * damp1D, float * damp2D, float * damp3D, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz, int nbzu, int nb);
__device__ float get_boundary_damper(float * damp1D, float * damp2D, float * damp3D, int i, int j, int k, int nxx, int nyy, int nzz, int nb, int nbzu);

# endif