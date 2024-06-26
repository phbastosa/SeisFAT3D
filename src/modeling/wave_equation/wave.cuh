# ifndef WAVE_CUH
# define WAVE_CUH

# include "../modeling.hpp"

class Wave : public Modeling
{
private:

    void set_outputs();
    void set_specifics();
    void define_cerjan_dampers();

protected:

    float dt, fmax, pabc;   

    int nt, nabc;
    int snap_index; 
    int total_snaps;

    float * damp1D = nullptr;
    float * damp2D = nullptr;
    float * damp3D = nullptr;

    float * wavelet = nullptr;

    int * grid_node_x = nullptr;
    int * grid_node_y = nullptr;
    int * grid_node_z = nullptr;

    float * snapshot = nullptr; 
    float * seismogram = nullptr;

    void display_progress();
    void define_common_wavelet();
    void define_staggered_wavelet();
    void define_grid_nodes_position();

    void get_seismogram();
    void get_receiver_output();
    void get_wavefield_output();

    virtual void set_models() = 0;
    virtual void set_volumes() = 0;
    
    virtual void initialization() = 0;

public: 

    virtual void set_forward_solver() = 0; 
    virtual void free_space() = 0;
};

__global__ void compute_seismogram(float * seismogram, float * P, int * rx, int * ry, int * rz, int total_nodes, int nxx, int nzz, int nt, int time_id);
__device__ float get_boundary_damper(float * damp1D, float * damp2D, float * damp3D, int i, int j, int k, int nxx, int nyy, int nzz, int nb);

# endif