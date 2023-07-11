# ifndef ACOUSTIC_HPP
# define ACOUSTIC_HPP

# include "../low_frequency.cuh" 

class Acoustic : public Low_Frequency
{
private:


protected:

    virtual void set_modeling_message() = 0;
    virtual void set_model_parameters() = 0;
    virtual void set_wavefields() = 0;
    
    void set_wavelet();

public:

    virtual void forward_solver() = 0;
    virtual void initial_setup() = 0;
    virtual void free_space() = 0;
};

# endif