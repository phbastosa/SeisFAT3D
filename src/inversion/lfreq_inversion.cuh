# ifndef LFREQ_INVERSION_CUH
# define LFREQ_INVERSION_CUH

# include "inversion.hpp"

# include "../modeling/lfreq_modeling.cuh"

class lfreq_Inversion : public Inversion
{
private:

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