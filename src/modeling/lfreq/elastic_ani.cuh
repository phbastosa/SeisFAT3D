# ifndef ELASTIC_VTI_CUH
# define ELASTIC_VTI_CUH

# include "elastic.cuh"

# include "../hfreq/eikonal_vti.cuh" 

class Elastic_VTI : public Elastic
{
private:

    float * T = nullptr;
    float * TT = nullptr;
    float * d_T = nullptr;

    float * E = nullptr;
    float * D = nullptr;
    float * G = nullptr;
    float * B = nullptr;

    float * C11 = nullptr;
    float * C33 = nullptr;
    float * C55 = nullptr;
    float * C66 = nullptr;
    float * C13 = nullptr;

    float * d_B = nullptr;
    float * d_C11 = nullptr;
    float * d_C33 = nullptr;
    float * d_C55 = nullptr;
    float * d_C66 = nullptr;
    float * d_C13 = nullptr;

    Modeling * eikonal = nullptr;

    void set_conditions();

public:

    void forward_solver();
};

__global__ void compute_pressure(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * P, float * T, 
                                 float * C11, float * C33, float * C55, float * C66, float * C13, float * wavelet, int sIdx, int sIdy, int sIdz, 
                                 int tId, int tlag, int nt, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz);

# endif