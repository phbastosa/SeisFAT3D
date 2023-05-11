# ifndef ACOUSTIC_HPP
# define ACOUSTIC_HPP

# include "../modeling.hpp"

class Acoustic : public Modeling
{

    
public:
    
    void initial_setup();
    void set_components();
    void forward_solver();
    void build_outputs();    
    void free_space();
};

# endif