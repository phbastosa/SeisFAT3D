# ifndef ACOUSTIC_CUH
# define ACOUSTIC_CUH

# include "../wave.hpp"

class Acoustic : public Wave
{
private:

    float * Vx = nullptr;
    float * Vy = nullptr;
    float * Vz = nullptr;

    void set_models();
    void set_volumes();
    void initialization();

protected:


public: 

    void set_forward_solver();

};

# endif