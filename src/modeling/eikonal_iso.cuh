# ifndef EIKONAL_ISO_CUH
# define EIKONAL_ISO_CUH

# include "eikonal.cuh"

class Eikonal_ISO : public Eikonal
{
private:

    void set_properties();
    void set_conditions();

public:

    void forward_solver();
};

# endif