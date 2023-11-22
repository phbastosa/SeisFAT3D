# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include <chrono>
# include <cuda_runtime.h>
# include <sys/resource.h>

# include "../geometry/regular/regular.hpp"
# include "../geometry/circular/circular.hpp"

class Eikonal
{
    int RAM, vRAM, ivRAM;

    void get_RAM_usage();
    void get_GPU_usage();
    void get_GPU_initMem();

    void check_geometry_overflow();

    std::chrono::_V2::system_clock::time_point ti, tf;

protected:

    int sidx, sidy, sidz;

    int receiver_output_samples;
    int wavefield_output_samples;

    bool export_receiver_output;
    bool export_wavefield_output;
    
    std::string eikonal_method;
    std::string eikonal_message;

    std::string receiver_output_file;
    std::string wavefield_output_file;
    std::string receiver_output_folder;
    std::string wavefield_output_folder;

    void set_outputs();
    void set_boundaries();
    void set_velocity_model();
    void set_slowness_model();
    void set_general_parameters();
    void set_acquisition_geometry();

    void get_travel_times();
    void get_first_arrivals();

    virtual void initialization() = 0;
    virtual void set_eikonal_volumes() = 0;
    virtual void set_specific_boundary() = 0;

public: 

    int shot_id;
    int total_shots;
    int total_nodes;
    
    std::string file;

    Geometry * geometry;

    float dx, dy, dz;
    int nbxl, nbxr, nbyl; 
    int nbyr, nbzu, nbzd;

    int blocksPerGrid;
    int threadsPerBlock;

    int source_id, time_id;
    int nx, ny, nz, nPoints;
    int nxx, nyy, nzz, volsize;

    float * V = nullptr;
    float * S = nullptr;
    float * T = nullptr;

    float * receiver_output = nullptr;
    float * wavefield_output = nullptr;

    void expand_boundary(float * input, float * output);
    void reduce_boundary(float * input, float * output);

    virtual void forward_solver() = 0;
    virtual void free_space() = 0;

    void set_parameters(); 
    void info_message();
    void initial_setup();
    void build_outputs();
    void export_outputs();

    void set_runtime();
    void get_runtime();
};

# endif