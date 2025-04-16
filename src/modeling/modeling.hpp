# ifndef MODELING_HPP
# define MODELING_HPP

# include <cuda_runtime.h>

# include "../geometry/geometry.hpp"

typedef unsigned short int uintc; 

class Modeling
{
private:

protected:
    
    float compress;
    float conversion;

    virtual void set_boundaries() = 0;
    virtual void set_specifications() = 0;

    void compression(float * input, uintc * output, int N, float &max_value, float &min_value, int levels);

public:

    float dx, dy, dz, dt;
    int nxx, nyy, nzz, volsize;
    int nt, nx, ny, nz, nb, nPoints;
    int srcId, recId, sIdx, sIdy, sIdz;

    float * S = nullptr;
    float * T = nullptr;
    float * Vp = nullptr;
    float * Vs = nullptr;
    float * Ro = nullptr;

    float * synthetic_data = nullptr;

    int max_spread;
    Geometry * geometry;
    
    std::string parameters;
    std::string data_folder;
    std::string modeling_type;
    std::string modeling_name;

    void set_parameters();
    void show_information();    

    void expand_boundary(float * input, float * output);
    void reduce_boundary(float * input, float * output);

    virtual void propagation() = 0;    
    virtual void initialization() = 0;
    virtual void forward_solver() = 0;

    virtual void export_synthetic_data() = 0;
};

# endif