# ifndef SHOOTING_HPP
# define SHOOTING_HPP


# include "../ray_tracing.hpp"

class Shooting : public Ray_Tracing
{
private:



protected:

    void set_modeling_message();
    void set_preconditioners();

public:

    void initial_setup();
    void forward_solver();
    void free_space();
};

# endif





