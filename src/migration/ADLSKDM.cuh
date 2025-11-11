# ifndef ADLSKDM_CUH
# define ADLSKDM_CUH

# include "LSKDM.cuh"

class ADLSKDM : public LSKDM
{
    void set_migration();
    void perform_forward();
    void perform_adjoint();

    void perform_adjoint_gradient();
    void perform_forward_direction();
};

# endif