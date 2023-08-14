# ifndef SCALAR_CUH
# define SCALAR_CUH

# include "../wave.cuh"

class Scalar : public Wave
{
protected:

    float * dtVp2 = nullptr;
    float * U_pre = nullptr;
    float * U_pas = nullptr;

    void set_wavelet();
    void set_model_boundaries();
    void set_modeling_volumes();
    
public:

    void set_parameters();
    void info_message();
    void forward_solver();
    void initial_setup();
    void free_space();
};

__global__ void compute_pressure(float * Pressure, float * U_pre, float * U_pas, float * dtVp2, float * damp1D, float * damp2D, float * damp3D, float * wavelet, int source_id, int time_id, float dx, float dy, float dz, int nxx, int nyy, int nzz, int nb, int nbzu);
__global__ void update_pressure(float * Pressure, float * U_pre, float * U_pas, int volsize);
__device__ float get_boundary_damper(float * damp1D, float * damp2D, float * damp3D, int i, int j, int k, int nxx, int nyy, int nzz, int nb, int nbzu);

# endif