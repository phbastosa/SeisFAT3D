# ifndef SCALAR_HPP
# define SCALAR_HPP

# include "../modeling.hpp"

class Scalar : public Modeling
{
protected:    

    int nxx, nyy, nzz, nb;    

    float * damp1D = nullptr;
    float * damp2D = nullptr;
    float * damp3D = nullptr;

    float * U_pre = nullptr;
    float * U_pas = nullptr;
    float * U_fut = nullptr;

    

public:
    
    void set_parameters();
    void set_components();
    void forward_solver();
    void initial_setup();
    void build_outputs();    
    void free_space();
};

# endif