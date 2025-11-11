# ifndef IDLSKDM_CUH
# define IDLSKDM_CUH

# include "LSKDM.cuh"

class IDLSKDM : public LSKDM
{
    void set_migration();
    void perform_forward();
    void perform_adjoint();

    void perform_adjoint_gradient();
    void perform_forward_direction();
};


# endif