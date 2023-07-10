# ifndef SCALAR_ISOTROPIC_HPP
# define SCALAR_ISOTROPIC_HPP

# include "../scalar.hpp"

class Scalar_Isotropic : public Scalar
{
private:


protected:

    float * Vp = nullptr;

    float * dtVp2 = nullptr;
    float * U_pre = nullptr;
    float * U_pas = nullptr;

    void set_modeling_message();
    void set_model_parameters();
    void set_wavefields();
    
public:

    void forward_solver();
    void initial_setup();
    void free_space();
};

__global__ void apply_wavelet(float * Pressure, float * wavelet, int sId, int time_id, float dx, float dy, float dz);
__global__ void compute_pressure(float * Pressure, float * U_pre, float * U_pas, float * dtVp2, float * damp1D, float * damp2D, float * damp3D, float dx, float dy, float dz, int nxx, int nyy, int nzz, int nb, int nbzu);
__global__ void update_pressure(float * Pressure, float * U_pre, float * U_pas, int volsize);

# endif