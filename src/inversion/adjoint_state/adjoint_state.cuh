# ifndef ADJOINT_STATE_CUH
# define ADJOINT_STATE_CUH

# include "../tomography.hpp"

class Adjoint_State : public Tomography
{
private:

    int totalLevels;
    int nSweeps, meshDim;

    float cell_volume;

    float * d_T = nullptr;

    float * source = nullptr;       
    float * adjoint = nullptr;
    
    float * d_source = nullptr;
    float * d_adjoint = nullptr;

    int iDivUp(int a, int b);

    void set_specific_parameters();
    void gradient_preconditioning();
    void apply_inversion_technique();

public:

    void optimization();
};

# endif

__global__ void adjoint_state_kernel(float * adjoint, float * source, float * T, int level, int xOffset, int yOffset, int xSweepOffset, int ySweepOffset, 
                                     int zSweepOffset, int nxx, int nyy, int nzz, float dx, float dy, float dz);


