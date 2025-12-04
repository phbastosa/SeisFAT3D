# ifndef EIKONAL_ISO_CUH
# define EIKONAL_ISO_CUH

# include "modeling.cuh"

class Eikonal_ISO : public Modeling
{
public:

    void set_conditions();
    void time_propagation();
};

# endif