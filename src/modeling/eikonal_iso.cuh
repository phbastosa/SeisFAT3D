# ifndef EIKONAL_ISO_CUH
# define EIKONAL_ISO_CUH

# include "modeling.cuh"

class Eikonal_ISO : public Modeling
{
private:

    void set_properties();
    void set_conditions();

public:

    void time_propagation();
};

# endif