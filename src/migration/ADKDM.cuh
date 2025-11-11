# ifndef ADKDM_CUH
# define ADKDM_CUH

# include "KDM.cuh"

class ADKDM : public KDM
{
    void set_migration();
    void perform_forward();
    void perform_adjoint();
};

# endif
