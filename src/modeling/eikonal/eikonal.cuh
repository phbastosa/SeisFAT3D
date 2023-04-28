# ifndef EIKONAL_CUH
# define EIKONAL_CUH

# include <cuda.h>
# include <cuda_runtime.h>

# include "../modeling.hpp"

# define INF 1e6f
# define EPS 1e-6f
# define BLOCK_LENGTH 4

# define MEM(index) _mem[index]
# define SOL(i,j,k) _sol[i][j][k]
# define SPD(i,j,k) _spd[i][j][k]

class Eikonal : public Modeling
{
private:

    float * vp = nullptr;

    int pdx, pdy, pdz; 
    int nxx, nyy, nzz;
    int nbx, nby, nbz;

    int nblk;
    int volsize;
    int blksize;
    int nit = 10;

    float * h_time;
    float * h_slow;

    uint * h_list;
    bool * h_mask;
    bool * h_listed;
    bool * h_listVol;

    void expand_model();

    void get_travelTimes();
    void get_firstArrivals();

public:

    float * T = nullptr;

    void initial_setup();
    void set_components();
    void forward_solver();
    void build_outputs();    

    void set_parameters(std::string file);
};

__global__ void run_solver(float* spd, bool* mask, const float *sol_in, float *sol_out, bool *con, uint* list, int xdim, int ydim, int zdim, float dh, int nIter, uint nActiveBlock);
__global__ void run_reduction(bool *con, bool *listVol, uint *list, uint nActiveBlock);
__global__ void run_check_neighbor(float* spd, bool* mask, const float *sol_in, float *sol_out, bool *con, uint* list, int xdim, int ydim, int zdim, float dh, uint nActiveBlock, uint nTotalBlock);
__device__ float get_time_eikonal(float a, float b, float c, float h, float s);

# endif
