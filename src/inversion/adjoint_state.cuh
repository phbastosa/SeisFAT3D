# ifndef ADJOINT_STATE_CUH
# define ADJOINT_STATE_CUH

# include "tomography.hpp"

class Adjoint_State : public Tomography
{
private:

    int total_levels;
    int nSweeps, meshDim;
    int nThreads, nBlocks;

    float cell_area;
    float aperture_x;
    float aperture_y;

    float * m = nullptr;
    float * v = nullptr;

    float * d_T = nullptr;

    float * d_source_grad = nullptr;
    float * d_source_comp = nullptr;

    float * d_adjoint_grad = nullptr;
    float * d_adjoint_comp = nullptr;

    float * h_source_grad = nullptr;
    float * h_source_comp = nullptr;

    float * h_adjoint_grad = nullptr;
    float * h_adjoint_comp = nullptr;

    float * gradient = nullptr;

    void initialization();
    void set_specifications();

    void apply_inversion_technique();

    int iDivUp(int a, int b);

public:

    void optimization();
};

__global__ void adjoint_state_kernel(float * T, float * adjoint_grad, float * adjoint_comp, float * source_grad, float * source_comp, int level, int xOffset, 
                                     int yOffset, int xSweepOffset, int ySweepOffset, int zSweepOffset, int nxx, int nyy, int nzz, float dx, float dy, float dz);

# endif
