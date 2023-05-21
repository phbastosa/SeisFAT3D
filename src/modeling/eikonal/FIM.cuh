# ifndef EIKONAL_FIM_CUH
# define EIKONAL_FIM_CUH

# include "eikonal.cuh"

# define INF 1e6f
# define EPS 1e-6f
# define BLOCK_LENGTH 4

# define MEM(index) _mem[index]
# define SOL(i,j,k) _sol[i][j][k]
# define SPD(i,j,k) _spd[i][j][k]

class Eikonal_fim : public Eikonal
{
    int nbx, nby, nbz;
    int pdx, pdy, pdz;
    int nblk, blksize;
    uint nActiveBlock;

    uint * h_list = nullptr;
    bool * h_mask = nullptr; 
    bool * h_listed = nullptr;
    bool * h_listVol = nullptr; 
    float * h_slow = nullptr; 
    float * h_time = nullptr; 
    
    uint * d_list = nullptr;
    bool * d_mask = nullptr;
    bool * d_con = nullptr;
    bool * d_listVol = nullptr;
    float * d_slow = nullptr;
    float * d_time = nullptr; 
    float * t_time = nullptr;

    void expansion();
    void reduction();
    void parameters();
    void components();

public:

    void free_space();
    void initial_setup();
    void forward_solver();
};

__global__ void run_solver(float* spd, bool* mask, const float *sol_in, float *sol_out, bool *con, uint* list, int xdim, int ydim, int zdim, float dh, int nIter, uint nActiveBlock);
__global__ void run_reduction(bool *con, bool *listVol, uint *list, uint nActiveBlock);
__global__ void run_check_neighbor(float* spd, bool* mask, const float *sol_in, float *sol_out, bool *con, uint* list, int xdim, int ydim, int zdim, float dh, uint nActiveBlock, uint nTotalBlock);
__device__ float get_time_eikonal(float a, float b, float c, float h, float s);

# endif