# ifndef WAVEFORM_HPP
# define WAVEFORM_HPP

# include "../inversion.hpp"

class Waveform : public Inversion
{
private:


protected:


public:

    void optimization();
    void model_update();
    void export_results();
    void set_parameters();
    void import_obs_data();
    void forward_modeling();
    void check_convergence();
};

# endif
