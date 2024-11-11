# ifndef ELASTIC_ISO_CUH
# define ELASTIC_ISO_CUH

# include "elastic.cuh"

class Elastic_Iso : public Elastic
{
private:

    float * M = nullptr;
    float * L = nullptr;
    float * B = nullptr;
    float * P = nullptr;

    float * d_M = nullptr;
    float * d_L = nullptr;
    float * d_B = nullptr;
    float * d_P = nullptr;

    float * d_Vx = nullptr;
    float * d_Vy = nullptr;
    float * d_Vz = nullptr;

    float * d_Txx = nullptr;
    float * d_Tyy = nullptr;
    float * d_Tzz = nullptr;
    
    float * d_Txz = nullptr;
    float * d_Tyz = nullptr;
    float * d_Txy = nullptr;

    void set_properties();
    void set_conditions();

public:

    void initialization();
    void forward_solver();
};

__global__ void compute_velocity(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * B,  float * damp1D, float * damp2D, float * damp3D, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz, int nb);
__global__ void compute_pressure(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * P, float * M, float * L, float * wavelet, int sIdx, int sIdy, int sIdz, int tId, int nt, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz);

__global__ void compute_seismogram(float * P, int * rIdx, int * rIdy, int * rIdz, float * seismogram, int spread, int tId, int tlag, int nt, int nxx, int nzz);

__device__ float get_boundary_damper(float * damp1D, float * damp2D, float * damp3D, int i, int j, int k, int nxx, int nyy, int nzz, int nabc);

# endif