# ifndef INVERSION_HPP
# define INVERSION_HPP

# include "../modeling/modeling.hpp"

class Inversion
{
protected:

    Modeling * modeling;

public:

    virtual void import_obs_data() = 0;
    virtual void forward_modeling() = 0;
    virtual void compute_residuals() = 0;
    virtual void check_convergence() = 0;
    virtual void optimization() = 0;
    virtual void model_update() = 0;

    virtual void set_parameters(std::string file) = 0;    

    void export_results();
};

# endif