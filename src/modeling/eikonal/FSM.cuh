# ifndef EIKONAL_FSM_CUH
# define EIKONAL_FSM_CUH

# include "eikonal.hpp"

class Eikonal_fsm : public Eikonal
{
private:

    float * d_T = nullptr;
    float * d_S = nullptr;

    int * d_sgnv = nullptr;
    int * d_sgnt = nullptr;

    int threadsPerBlock;
    int nSweeps, meshDim;
    int totalLevels, padb;
    float dh2i, dh4i, dsum;

    void expansion();
    void reduction();
    void parameters();
    void components();

    int iDivUp(int a, int b); 

public:

    void free_space();
    void initial_setup();
    void forward_solver();
};

__global__ void fast_sweeping_kernel(float * S, float * T, int * sgnt, int * sgnv, int sc, int nSweeps, 
                                     int level, int xOffset, int yOffset, int xSweepOffset, int ySweepOffset, 
                                     int zSweepOffset, int nxx, int nyy, int nzz, float dh, float dh2i, 
                                     float dh4i, float dsum);

# endif