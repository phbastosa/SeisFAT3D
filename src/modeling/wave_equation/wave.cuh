# ifndef LOW_FREQUENCY_HPP
# define LOW_FREQUENCY_HPP

# include "../modeling.hpp" 

class Wave : public Modeling
{
private:


protected:

    int time_id, nt, sId;
    int nsnap, dsnap, isnap;

    float fmax, tlag, amp;
    float dt, tmax, pi, fc;

    bool free_surface;
    bool import_wavelet;
    std::string wavelet_file;

    float * damp1D = nullptr;
    float * damp2D = nullptr;
    float * damp3D = nullptr;

    float * wavelet = nullptr;

    float * Pressure = nullptr;
    float * snapshot = nullptr;
    float * seismogram = nullptr;

    int * grid_node_x = nullptr;
    int * grid_node_y = nullptr;
    int * grid_node_z = nullptr;

    void set_specifications();
    void set_model_boundaries();
    void set_gridded_geometry();
    
    void set_dampers();
    void set_outputs();

    void show_progress();
    void get_snapshots();
    void get_seismogram();

    void set_modeling_message();

    virtual void set_model_parameters() = 0;
    virtual void set_wavefields() = 0;
    virtual void set_wavelet() = 0;

public:

    void build_outputs();

    virtual void forward_solver() = 0;
    virtual void initial_setup() = 0;
    virtual void free_space() = 0;
};

__global__ void compute_seismogram(float * Pressure, float * seismogram, int * rx, int * ry, int * rz, int total_nodes, int time_id, int nt, int nxx, int nzz);
__device__ float get_boundary_damper(float * damp1D, float * damp2D, float * damp3D, int i, int j, int k, int nxx, int nyy, int nzz, int nb, int nbzu);

# endif