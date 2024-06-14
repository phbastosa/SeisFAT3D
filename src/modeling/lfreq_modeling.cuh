# ifndef LFREQ_MODELING_CUH
# define LFREQ_MODELING_CUH

# include <cuda_runtime.h>

# include "modeling.hpp"

class lfreq_Modeling : public Modeling
{
private:

    int nt, nabc;
    int snap_index; 
    int total_snaps;

    float dt, fmax;   

    float * d_Vp = nullptr;
    float * d_Upas = nullptr;
    float * d_Upre = nullptr;
    float * d_Ufut = nullptr;

    float * damp1D = nullptr;
    float * damp2D = nullptr;
    float * damp3D = nullptr;

    float * wavelet = nullptr;

    int * grid_node_x = nullptr;
    int * grid_node_y = nullptr;
    int * grid_node_z = nullptr;

    float * seismogram = nullptr;

    void display_progress();
    void define_cerjan_dampers();
    void define_wavelet_signature();
    void define_grid_nodes_position();

    void set_outputs();

    void set_volumes();
    void set_specifics();
    void initialization();

    void get_receiver_output();
    void get_wavefield_output();
    
    void get_seismogram();

public: 

    void forward_propagation(); 
    void free_space();
};

__global__ void update_pressure(float * P, float * Pold, float * Pnew, int volsize);
__global__ void compute_pressure(float * P, float * U_old, float * U_new, float * V, float * damp1D, float * damp2D, float * damp3D, float * wavelet, int sId, int tId, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz, int nabc);
__global__ void compute_seismogram(float * seismogram, float * P, int * rx, int * ry, int * rz, int total_nodes, int nxx, int nzz, int nt, int time_id);
__device__ float get_boundary_damper(float * damp1D, float * damp2D, float * damp3D, int i, int j, int k, int nxx, int nyy, int nzz, int nb);

# endif
