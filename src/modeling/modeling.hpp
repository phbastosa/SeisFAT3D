# ifndef MODELING_HPP
# define MODELING_HPP

# include <string>
# include <vector> 
# include <iostream>

# include "geometry/geometry.hpp"
# include "geometry/regular/regular.hpp"
# include "geometry/circular/circular.hpp"

class Modeling
{
private:



protected:

    int receiver_output_samples;
    int wavefield_output_samples;

    float * receiver_output = nullptr;
    float * wavefield_output = nullptr;

    float dx, dy, dz;
    int nxx, nyy, nzz;
    int nx, ny, nz, nb;

    float * S = nullptr;
    float * T = nullptr;

public: 

    virtual void set_parameters(std::string file); 

    // virtual void set_components() = 0;

    // virtual void initial_setup() = 0;
    // virtual void forward_solver() = 0;

    void export_outputs();
};

# endif