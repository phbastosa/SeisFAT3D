# ifndef MODELING_HPP
# define MODELING_HPP

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

    std::chrono::_V2::system_clock::time_point ti, tf;

    void check_geometry_overflow();

protected:

    int receiver_output_samples;
    int wavefield_output_samples;

    bool export_receiver_output;
    bool export_wavefield_output;

    std::string receiver_output_file;
    std::string wavefield_output_file;
    std::string receiver_output_folder;
    std::string wavefield_output_folder;

    Geometry * geometry;

public: 

    int shot_id;

    int total_shots;
    int total_nodes;

    int nPoints;
    float dx, dy, dz;
    int nx, ny, nz, nt; 

    float * S = nullptr;
    float * V = nullptr;
    float * T = nullptr;
    
    float * receiver_output = nullptr;
    float * wavefield_output = nullptr;
    
    std::string file;

    void set_runtime();
    void get_runtime();
    void info_message();
    void export_outputs();

    virtual void set_parameters(); 
    virtual void set_components() = 0;
    virtual void forward_solver() = 0;
    virtual void initial_setup() = 0;
    virtual void build_outputs() = 0;
    virtual void free_space() = 0;
};

# endif