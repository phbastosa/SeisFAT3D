# ifndef ACOUSTIC_ISOTROPIC_HPP
# define ACOUSTIC_ISOTROPIC_HPP

# include "../acoustic.hpp" 

class Acoustic_Isotropic : public Acoustic
{
private:


protected:

    float * Vp = nullptr;
    float * Rho = nullptr;

    float * K = nullptr;
    float * B = nullptr;
    float * Vx = nullptr; 
    float * Vy = nullptr; 
    float * Vz = nullptr; 

    void set_modeling_message();
    void set_model_parameters();
    void set_wavefields();

public:

    void forward_solver();
    void initial_setup();
    void free_space();
};

# endif