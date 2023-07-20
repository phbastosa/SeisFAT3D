# ifndef ACOUSTIC_ISOTROPIC_HPP
# define ACOUSTIC_ISOTROPIC_HPP

# include "../acoustic.hpp" 

class Acoustic_Isotropic : public Acoustic
{
private:


protected:

    float * Vp = nullptr;
    float * Rho = nullptr;

    float * K = nullptr;
    float * B = nullptr;
    float * Vx = nullptr; 
    float * Vy = nullptr; 
    float * Vz = nullptr; 

    void set_model_parameters();
    void set_wavefields();

public:

    void forward_solver();
    void initial_setup();
    void free_space();
};

__global__ void compute_velocity(float * Pressure, float * Vx, float * Vy, float * Vz, float * B, float * wavelet, int sId, int time_id, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz);
__global__ void compute_pressure(float * Pressure, float * Vx, float * Vy, float * Vz, float * K, float * damp1D, float * damp2D, float * damp3D, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz, int nbzu, int nb);
__device__ float get_boundary_damper(float * damp1D, float * damp2D, float * damp3D, int i, int j, int k, int nxx, int nyy, int nzz, int nb, int nbzu);

# endif