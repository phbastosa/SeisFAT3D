# ifndef EIKONAL_PAL_CUH
# define EIKONAL_PAL_CUH

# include "eikonal.hpp"

class Eikonal_pal : public Eikonal
{
    int padb;

    float * K = nullptr;
    
    float * d_K = nullptr;
    float * d_nK = nullptr;
    float * d_nT = nullptr;
    float * d_T = nullptr;
    float * d_S = nullptr;

    void expansion();
    void reduction();
    void parameters();
    void components();

public:

    void free_space();
    void initial_setup();
    void forward_solver();    
};

# endif

__global__ void equations(float * S, float * T, float * K, float * nT, float h, int nxx, int nyy, int nzz);
__global__ void wavefront(float * K, float * nK, int nxx, int nyy, int nzz);
__global__ void update(float * T, float * nT, float * K, float * nK, int N);
