# ifndef INVERSION_HPP
# define INVERSION_HPP

# include "../modeling/modeling.hpp"

class Inversion
{
protected:

    int n_data;
    int n_model;

    float tolerance;

    int iteration;
    int max_iteration;

    float * dobs = nullptr;
    float * dcal = nullptr;

    float * gradient = nullptr;
    
    std::string obs_data_folder;
    std::string obs_data_prefix;

    std::vector<float> residuo;

    Modeling * modeling;

public:

    std::string file;

    bool converged;

    virtual void set_parameters() = 0;    
    virtual void import_obs_data() = 0;
    virtual void forward_modeling() = 0;
    virtual void compute_residuals() = 0;
    virtual void check_convergence() = 0;
    virtual void optimization() = 0;
    virtual void model_update() = 0;

    void export_results();
};

# endif