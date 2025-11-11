# ifndef IDKDM_CUH
# define IDKDM_CUH

# include "KDM.cuh"

class IDKDM : public KDM
{
    void set_migration();

    void perform_forward();
    void perform_adjoint();
};

# endif
