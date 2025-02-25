# ifndef ELASTIC_ISO_CUH
# define ELASTIC_ISO_CUH

# include "elastic.cuh"

# include "../hfreq/eikonal_iso.cuh" 

class Elastic_ISO : public Elastic
{
private:

    float * M = nullptr;
    float * L = nullptr;
    float * B = nullptr;

    float * d_M = nullptr;
    float * d_L = nullptr;
    float * d_B = nullptr;

    void set_conditions();

public:

    void propagation();
};

__global__ void compute_velocity_ssg(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * B, float * T,  
                                     float * damp1D, float * damp2D, float * damp3D, float * wavelet, float dx, float dy, float dz, float dt, int tId, int tlag, int sIdx, 
                                     int sIdy, int sIdz, int nxx, int nyy, int nzz, int nb, int nt);

__global__ void compute_pressure_ssg(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * P, float * M, 
                                     float * L, float * T, int tId, int tlag, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz);

# endif