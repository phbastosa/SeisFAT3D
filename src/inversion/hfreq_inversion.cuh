# ifndef HFREQ_INVERSION_CUH
# define HFREQ_INVERSION_CUH

# include "inversion.hpp"

# include "../modeling/hfreq_modeling.cuh"

class hfreq_Inversion : public Inversion
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

    void get_objective_function();

    void set_forward_modeling();
    void set_inversion_volumes();

    void extract_calculated_data();
    void adjoint_propagation();

    void update_specifications();

public:

    void import_obs_data();

};

# endif

__global__ void adjoint_state_kernel(float * adjoint, float * source, float * T, int level, int xOffset, int yOffset, int xSweepOffset, int ySweepOffset, 
                                     int zSweepOffset, int nxx, int nyy, int nzz, float dx, float dy, float dz);
