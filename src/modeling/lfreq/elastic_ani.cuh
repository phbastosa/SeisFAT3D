# ifndef ELASTIC_ANI_CUH
# define ELASTIC_ANI_CUH

# include "elastic.cuh"

# include "../hfreq/eikonal_ani.cuh" 

class Elastic_ANI : public Elastic
{
private:

    float * T = nullptr;
    float * TT = nullptr;
    float * d_T = nullptr;

    float * B = nullptr;

    float * C11 = nullptr;
    float * C12 = nullptr;
    float * C13 = nullptr;
    float * C14 = nullptr;
    float * C15 = nullptr;
    float * C16 = nullptr;

    float * C22 = nullptr;
    float * C23 = nullptr;
    float * C24 = nullptr;
    float * C25 = nullptr;
    float * C26 = nullptr;

    float * C33 = nullptr;
    float * C34 = nullptr;
    float * C35 = nullptr;
    float * C36 = nullptr;

    float * C44 = nullptr;
    float * C45 = nullptr;
    float * C46 = nullptr;

    float * C55 = nullptr;
    float * C56 = nullptr;

    float * C66 = nullptr;

    float * d_B = nullptr;

    float * d_C11 = nullptr;
    float * d_C12 = nullptr;
    float * d_C13 = nullptr;
    float * d_C14 = nullptr;
    float * d_C15 = nullptr;
    float * d_C16 = nullptr;

    float * d_C22 = nullptr;
    float * d_C23 = nullptr;
    float * d_C24 = nullptr;
    float * d_C25 = nullptr;
    float * d_C26 = nullptr;

    float * d_C33 = nullptr;
    float * d_C34 = nullptr;
    float * d_C35 = nullptr;
    float * d_C36 = nullptr;

    float * d_C44 = nullptr;
    float * d_C45 = nullptr;
    float * d_C46 = nullptr;

    float * d_C55 = nullptr;
    float * d_C56 = nullptr;

    float * d_C66 = nullptr;

    void set_conditions();

public:

    void propagation();
};

__global__ void compute_pressure(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * P, float * T, 
                                 float * C11, float * C12, float * C13, float * C14, float * C15, float * C16, float * C22, float * C23, float * C24, float * C25, float * C26,  
                                 float * C33, float * C34, float * C35, float * C36, float * C44, float * C45, float * C46, float * C55, float * C56, float * C66, float * wavelet, 
                                 int sIdx, int sIdy, int sIdz, int tId, int tlag, int nt, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz);

# endif