# ifndef ELASTIC_ISO_CUH
# define ELASTIC_ISO_CUH

# include "elastic.cuh"

# include "../hfreq/eikonal_iso.cuh" 

class Elastic_ISO : public Elastic
{
private:

    float * T = nullptr;
    float * TT = nullptr;
    float * d_T = nullptr;

    float * M = nullptr;
    float * L = nullptr;
    float * B = nullptr;

    float * d_M = nullptr;
    float * d_L = nullptr;
    float * d_B = nullptr;

    Modeling * eikonal = nullptr;

    void set_conditions();

public:

    void propagation();
    void forward_solver();
};

__global__ void compute_pressure(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * P, float * M, float * L, float * T, float * wavelet, int sIdx, int sIdy, int sIdz, int tId, int tlag, int nt, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz);

# endif