# ifndef SCALAR_CUH
# define SCALAR_CUH

# include "../wave.hpp"

class Scalar : public Wave
{
private:

    float * Pold = nullptr;
    float * Pnew = nullptr;

    void set_models();
    void set_volumes();
    void initialization();

protected:


public: 

    void set_forward_solver();
    
};

# endif