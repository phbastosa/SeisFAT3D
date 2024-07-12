# ifndef KIRCHHOFF_CUH
# define KIRCHHOFF_CUH

# include "../migration.hpp"

# include "../../modeling/eikonal_equation/isotropic/ultimate_FSM.cuh"

class Kirchhoff : public Migration
{
private:

    void set_modeling_type();

protected:


     
public:
    
    void image_building();
};

# endif