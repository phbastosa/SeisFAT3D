# ifndef SCALAR_HPP
# define SCALAR_HPP

# include "../modeling.hpp"

class Scalar : public Modeling
{


public:
    
    void initial_setup();
    void set_components();
    void forward_solver();
    void build_outputs();    
    void free_space();
};

# endif