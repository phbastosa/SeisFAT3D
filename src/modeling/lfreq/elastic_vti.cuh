# ifndef ELASTIC_VTI_CUH
# define ELASTIC_VTI_CUH

# include "elastic.cuh"

class Elastic_VTI : public Elastic
{
private:

    float * E = nullptr;
    float * D = nullptr;
    float * G = nullptr;

    float * d_B = nullptr;
    float * d_C11 = nullptr;
    float * d_C33 = nullptr;
    float * d_C55 = nullptr;
    float * d_C66 = nullptr;
    float * d_C13 = nullptr;

    void set_conditions();

public:

    void forward_solver();
};

__global__ void compute_pressure(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, 
                                 float * P, float * C11, float * C33, float * C55, float * C66, float * C13, float * wavelet, int sIdx, int sIdy, int sIdz, 
                                 int tId, int nt, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz);

# endif