# ifndef WAVE_HPP
# define WAVE_HPP

# include "../modeling.hpp"

class Wave : public Modeling
{
private:

    void set_specifics();

    void define_cerjan_dampers();

    void define_seismogram();
    void define_snapshots();

protected:

    float dt, fmax, pabc;   

    int nt, nabc, total_snaps;

    int time_index;

    float * damp1D = nullptr;
    float * damp2D = nullptr;
    float * damp3D = nullptr;

    float * wavelet = nullptr;

    float * snapshots = nullptr;
    float * seismogram = nullptr;

    void define_common_wavelet();
    void define_staggered_wavelet();

    virtual void set_models() = 0;
    virtual void set_volumes() = 0;
    
    virtual void initialization() = 0;

public: 

    void set_forward_solver() = 0; 

};

# endif