# ifndef WAVEFORM_HPP
# define WAVEFORM_HPP

# include "../inversion.hpp"

class Waveform : public Inversion
{

public:
    void import_obs_data();
    void forward_modeling();
    void compute_residuals();
    void check_convergence();
    void optimization();
    void model_update();

    void set_parameters(std::string file);
};

# endif