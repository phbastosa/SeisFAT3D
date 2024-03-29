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
    void free_space();
};

__global__ void compute_velocity(float * P, float * Vx, float * Vy, float * Vz, float * B, float * wavelet, int sId, int tId, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz);
__global__ void compute_pressure(float * P, float * Vx, float * Vy, float * Vz, float * K, float * damp1D, float * damp2D, float * damp3D, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz, int nabc);


# endif