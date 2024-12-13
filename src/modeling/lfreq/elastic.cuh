# ifndef ELASTIC_CUH
# define ELASTIC_CUH

# include "../modeling.hpp"

class Elastic : public Modeling
{
private:

    void set_wavelet();
    void set_boundaries();
    void set_properties();
    void set_specifications();

protected:

    float fmax, bd;

    int tlag, nThreads;
    int sBlocks, nBlocks;

    float * P = nullptr;
    float * d_P = nullptr;

    float * d_Vx = nullptr;
    float * d_Vy = nullptr;
    float * d_Vz = nullptr;

    float * d_Txx = nullptr;
    float * d_Tyy = nullptr;
    float * d_Tzz = nullptr;
    
    float * d_Txz = nullptr;
    float * d_Tyz = nullptr;
    float * d_Txy = nullptr;

    float * d1D = nullptr;
    float * d2D = nullptr;
    float * d3D = nullptr;

    int * rIdx = nullptr;
    int * rIdy = nullptr;
    int * rIdz = nullptr;

    int * current_xrec = nullptr;
    int * current_yrec = nullptr;
    int * current_zrec = nullptr;

    float * wavelet = nullptr;

    float * seismogram = nullptr;

    virtual void set_conditions() = 0;

public:

    void initialization();
    
    virtual void forward_solver() = 0;

    void export_synthetic_data();
};

# endif