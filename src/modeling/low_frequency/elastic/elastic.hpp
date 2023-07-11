# ifndef ELASTIC_HPP
# define ELASTIC_HPP

# include "../low_frequency.cuh"

class Elastic : public Low_Frequency
{
private:

protected:

    virtual void set_model_parameters() = 0;
    virtual void set_wavefields() = 0;
    
    void set_wavelet();

public:

    virtual void initial_setup() = 0;
    virtual void forward_solver() = 0;
    virtual void free_space() = 0;
};


# endif