# ifndef INVERSION_HPP
# define INVERSION_HPP

# include "../utils/input_output/io.hpp"

class Inversion
{
protected:

    float * dobs = nullptr;
    float * dcal = nullptr;

    int n_data, n_model;
    int max_iteration;
    int iteration;

    float tolerance;

    std::vector<float> residuo;

    std::string obs_data_folder;
    std::string obs_data_prefix;

    std::string convergence_map_folder;
    std::string estimated_model_folder;

    std::string inversion_method;

    void set_general_parameters();

public:

    bool converged;

    std::string file;

    virtual void set_parameters() = 0;
    virtual void import_obs_data() = 0;

    virtual void forward_modeling() = 0;
    virtual void check_convergence() = 0;

    virtual void optimization() = 0;
    virtual void model_update() = 0;
    virtual void export_results() = 0;
};

# endif
