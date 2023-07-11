# ifndef MODELING_HPP
# define MODELING_HPP

# include <omp.h>
# include <chrono>
# include <cuda_runtime.h>
# include <sys/resource.h>

# include "../geometry/geometry.hpp"
# include "../geometry/regular/regular.hpp"
# include "../geometry/circular/circular.hpp"

class Modeling
{
private:

    int RAM, vRAM, ivRAM;

    void get_RAM_usage();
    void get_GPU_usage();
    void get_GPU_initMem();

    void check_geometry_overflow();

    std::chrono::_V2::system_clock::time_point ti, tf;

protected:

    int blocksPerGrid;
    int threadsPerBlock;

    int nx, ny, nz, nPoints;
    int nxx, nyy, nzz, nb, volsize;
    int nbxl, nbxr, nbyl, nbyr, nbzu, nbzd;

    float dx, dy, dz;

    int receiver_output_samples;
    int wavefield_output_samples;

    bool export_receiver_output;
    bool export_wavefield_output;

    float * receiver_output = nullptr;
    float * wavefield_output = nullptr;
    
    std::string modeling_method;

    std::string receiver_output_file;
    std::string wavefield_output_file;
    std::string receiver_output_folder;
    std::string wavefield_output_folder;

    Geometry * geometry;

    virtual void set_specifications() = 0;
    virtual void set_modeling_message() = 0;

    void expand_boundary(float * input, float * output);

public: 

    int shot_id;
    int time_id;
    int total_shots;
    int total_nodes;
    
    std::string file;

    void set_runtime();
    void get_runtime();
    void info_message();
    void set_parameters(); 
    void export_outputs();

    virtual void forward_solver() = 0;
    virtual void initial_setup() = 0;
    virtual void build_outputs() = 0;
    virtual void free_space() = 0;
};

# endif