# ifndef ELASTIC_CUH
# define ELASTIC_CUH

# include "../wave.hpp"

class Elastic : public Wave
{
private:

    float * Vx = nullptr;
    float * Vy = nullptr;
    float * Vz = nullptr;

    float * Txx = nullptr;
    float * Tyy = nullptr;
    float * Tzz = nullptr;
    float * Txy = nullptr;
    float * Txz = nullptr;
    float * Tyz = nullptr;

    void set_models();
    void set_volumes();
    void initialization();

protected:


public: 

    void set_forward_solver();

};

# endif