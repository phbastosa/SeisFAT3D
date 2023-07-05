# ifndef ELASTIC_HPP
# define ELASTIC_HPP

# include "../modeling.hpp"

class Elastic : public Modeling
{
private:

    int nsnap, dsnap, isnap;
    int totalBlocks, sId;

    float * d_damp1D = nullptr;
    float * d_damp2D = nullptr;
    float * d_damp3D = nullptr;

    float * d_wavelet = nullptr;

    void set_model_parameters();
    void set_model_boundaries();
    void set_gridded_geometry();
    void set_wavefields();
    void show_progress();
    void set_dampers();
    void set_wavelet();
    void set_outputs();

    void specific_modeling_parameters();

protected:

    float * Vp = nullptr;
    float * Vs = nullptr;
    float * Rho = nullptr;

    float * h_B = nullptr;
    float * h_M = nullptr;
    float * h_L = nullptr;

    float * B = nullptr;
    float * L = nullptr;
    float * M = nullptr;

    float * Vx = nullptr;
    float * Vy = nullptr;
    float * Vz = nullptr;

    float * Txx = nullptr;
    float * Tyy = nullptr;
    float * Tzz = nullptr;
    float * Txz = nullptr;
    float * Tyz = nullptr;
    float * Txy = nullptr;

    float * Pressure = nullptr;
    float * snapshot = nullptr;
    float * seismogram = nullptr;

    int * grid_node_x = nullptr;
    int * grid_node_y = nullptr;
    int * grid_node_z = nullptr;

public:

    void initial_setup();
    void forward_solver();
    void build_outputs();
    void free_space();
};

__global__ void apply_wavelet(float * Txx, float * Tyy, float * Tzz, float * wavelet, int sId, int time_id, float dx, float dy, float dz);
__global__ void compute_velocity(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * B, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz);
__global__ void compute_stress(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * Pressure, float * M, float * L, float * damp1D, float * damp2D, float * damp3D, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz, int nb);
__global__ void get_seismogram(float * Pressure, float * seismogram, int * rx, int * ry, int * rz, int total_nodes, int time_id, int nt, int nxx, int nzz);

# endif