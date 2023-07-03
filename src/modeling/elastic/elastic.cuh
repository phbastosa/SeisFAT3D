# ifndef ELASTIC_HPP
# define ELASTIC_HPP

# include "../modeling.hpp"

class Elastic : public Modeling
{
private:

    void set_model_parameters();
    void set_model_boundaries();
    void set_wavefields();
    void set_outputs();

protected:


public:

    void initial_setup();
    void forward_solver();
    void build_outputs();
    void free_space();
};

# endif