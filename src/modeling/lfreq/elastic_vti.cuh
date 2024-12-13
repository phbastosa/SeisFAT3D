# ifndef ELASTIC_VTI_CUH
# define ELASTIC_VTI_CUH

# include "elastic.cuh"

class Elastic_VTI : public Elastic
{
private:

    float * E = nullptr;
    float * D = nullptr;
    float * G = nullptr;

    float * d_E = nullptr;
    float * d_D = nullptr;
    float * d_G = nullptr;

    void set_conditions();

public:

    void forward_solver();
};

// __global__ void compute_velocity(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * B,  float * damp1D, float * damp2D, float * damp3D, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz, int nb);
// __global__ void compute_pressure(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * P, float * M, float * L, float * wavelet, int sIdx, int sIdy, int sIdz, int tId, int nt, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz);

// __global__ void compute_seismogram(float * P, int * rIdx, int * rIdy, int * rIdz, float * seismogram, int spread, int tId, int tlag, int nt, int nxx, int nzz);

// __device__ float get_boundary_damper(float * damp1D, float * damp2D, float * damp3D, int i, int j, int k, int nxx, int nyy, int nzz, int nabc);

# endif