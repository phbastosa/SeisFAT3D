# ifndef ADJOINT_STATE_CUH
# define ADJOINT_STATE_CUH

# include "../tomography.hpp"

class Adjoint_State : public Tomography
{
private:

    int totalLevels;
    int nSweeps, meshDim;

    float * d_T = nullptr;

    int * d_sgnv = nullptr;
    int * d_sgnt = nullptr;

    float * source = nullptr;       
    float * adjoint = nullptr;
    
    float * d_source = nullptr;
    float * d_adjoint = nullptr;

    void initial_setup();
    void adjoint_state_solver();

    int iDivUp(int a, int b);

public:

    void optimization();
    void set_parameters();
    void forward_modeling();
};

# endif

__global__ void adjoint_state_kernel(float * adjoint, float * source, float * T, int level, int xOffset, int yOffset, 
                                     int xSweepOffset, int ySweepOffset, int zSweepOffset, int nxx, int nyy, int nzz, 
                                     float dx, float dy, float dz);


