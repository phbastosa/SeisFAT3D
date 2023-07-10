# ifndef ELASTIC_ISOTROPIC_HPP
# define ELASTIC_ISOTROPIC_HPP

# include "../elastic.hpp"

class Elastic_Isotropic : public Elastic
{
private:

    float * Vp;
    float * Vs;
    float * Rho;

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

protected:

    void set_modeling_message();
    void set_model_parameters();
    void set_wavefields();

public:

    void initial_setup();
    void forward_solver();
    void free_space();
};

__global__ void apply_wavelet(float * Txx, float * Tyy, float * Tzz, float * wavelet, int sId, int time_id, float dx, float dy, float dz);
__global__ void compute_velocity(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * B, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz);
__global__ void compute_stress(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * Pressure, float * M, float * L, float * damp1D, float * damp2D, float * damp3D, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz, int nb, int nbzu);

# endif