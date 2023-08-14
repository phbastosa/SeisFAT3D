# ifndef FAST_ITERATIVE_METHOD_HPP
# define FAST_ITERATIVE_METHOD_HPP

# include "../eikonal.hpp"

# include <cassert>

# define INF 1e6f
# define EPS 1e-6f
# define BLOCK_LENGTH 4

# define MEM(index) _mem[index]
# define SOL(i,j,k) _sol[i][j][k]
# define SPD(i,j,k) _spd[i][j][k]

class Fast_Iterative_Method : public Eikonal
{
private:

    int blksize;
    int nblk, nit;
    int nbx, nby, nbz;
    uint nActiveBlock;

    uint * h_list = nullptr;
    bool * h_listed = nullptr;
    bool * h_listVol = nullptr; 

    bool * h_mask = nullptr; 
    float * h_slow = nullptr; 
    float * h_time = nullptr; 
    
    bool * d_con = nullptr;
    uint * d_list = nullptr;
    bool * d_listVol = nullptr;

    bool * d_mask = nullptr;
    float * d_slow = nullptr;
    float * d_time = nullptr; 
    float * t_time = nullptr;

    void set_boundaries();
    void set_modeling_volumes();
    void check_spatial_spacing();

public:

    void set_parameters(); 

    void info_message();
    void initial_setup();
    void forward_solver();

    void free_space();
};

__global__ void run_solver(float* spd, bool* mask, const float *sol_in, float *sol_out, bool *con, uint* list, int xdim, int ydim, int zdim, float dh, int nIter, uint nActiveBlock);
__global__ void run_reduction(bool *con, bool *listVol, uint *list, uint nActiveBlock);
__global__ void run_check_neighbor(float* spd, bool* mask, const float *sol_in, float *sol_out, bool *con, uint* list, int xdim, int ydim, int zdim, float dh, uint nActiveBlock, uint nTotalBlock);
__device__ float get_time_eikonal(float a, float b, float c, float h, float s);

# endif
