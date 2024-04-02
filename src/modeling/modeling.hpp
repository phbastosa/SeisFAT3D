# ifndef MODELING_HPP
# define MODELING_HPP

# include <chrono>
# include <cuda_runtime.h>
# include <sys/resource.h>

# include "../geometry/regular/regular.hpp"
# include "../geometry/circular/circular.hpp"

class Modeling
{
private:

    std::chrono::system_clock::time_point ti, tf;

    int RAM, vRAM, ivRAM;

    void get_RAM_usage();
    void get_GPU_usage();
    void get_GPU_initMem();

    void set_generals();
    void set_geometry();

    void check_geometry_overflow();
    
protected:

    std::string type_name;
    std::string type_message;

    bool export_receiver_output;
    bool export_wavefield_output;

    int receiver_output_samples;
    int wavefield_output_samples;

    std::string receiver_output_file;
    std::string wavefield_output_file;
    std::string receiver_output_folder;
    std::string wavefield_output_folder;

    int sidx, sidy, sidz;

    void set_boundary();

    virtual void set_models() = 0;
    virtual void set_volumes() = 0;
    virtual void set_outputs() = 0;
    virtual void set_specifics() = 0;

    virtual void initialization() = 0;

    virtual void get_receiver_output() = 0;
    virtual void get_wavefield_output() = 0;

public:

    std::string file;

    int shot_index;
    int time_index;
    int source_index;    

    int total_shots;
    int total_nodes;

    int blocksPerGrid;
    int threadsPerBlock;

    float dx, dy, dz;
    int nx, ny, nz, nPoints;
    int nxx, nyy, nzz, volsize;

    int nbxl, nbxr, nbyl; 
    int nbyr, nbzu, nbzd;

    Geometry * geometry;

    float * S = nullptr;
    float * V = nullptr;
    float * K = nullptr;
    float * B = nullptr;
    float * M = nullptr;
    float * L = nullptr;

    float * T = nullptr;
    float * P = nullptr;

    float * receiver_output = nullptr;
    float * wavefield_output = nullptr;

    void expand_boundary(float * input, float * output);
    void reduce_boundary(float * input, float * output);

    void set_runtime();
    void get_runtime();

    void set_parameters(); 
    void get_information();
    void set_configuration();
    
    void export_outputs();

    virtual void set_forward_solver() = 0;
    virtual void free_space() = 0;    
};

# endif