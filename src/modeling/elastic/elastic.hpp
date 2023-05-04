# ifndef ELASTIC_HPP
# define ELASTIC_HPP

# include "../modeling.hpp"

class Elastic : public Modeling
{

public:
    
    void initial_setup();
    void set_components();
    void forward_solver();
    void build_outputs();    

};

# endif