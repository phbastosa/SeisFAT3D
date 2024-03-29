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
    void free_space();
};

__global__ void compute_pressure(float * P, float * U_old, float * U_new, float * V, float * damp1D, float * damp2D, float * damp3D, float * wavelet, int sId, int tId, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz, int nabc);
__global__ void update_pressure(float * P, float * Pold, float * Pnew, int volsize);

# endif