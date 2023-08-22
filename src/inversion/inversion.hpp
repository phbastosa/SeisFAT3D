# ifndef INVERSION_HPP
# define INVERSION_HPP

# include "../utils/input_output/io.hpp"

class Inversion
{
protected:

    int iteration;
    int max_iteration;
    int n_data, n_model;

    float tolerance;

    bool smooth;
    int samples;
    float stdv;

    float * dobs = nullptr;
    float * dcal = nullptr;

    float * dm = nullptr;
    float * model = nullptr;
    float * gradient = nullptr;

    bool write_model_per_iteration;
    bool write_gradient_per_iteration;

    std::vector<float> residuo;

    std::string obs_data_folder;
    std::string obs_data_prefix;

    std::string gradient_folder;
    std::string convergence_map_folder;
    std::string estimated_model_folder;

    std::string inversion_method;

    void set_general_parameters();
    void smoothing(float * input, float * output, int nx, int ny, int nz);

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
