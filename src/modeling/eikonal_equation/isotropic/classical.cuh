# ifndef CLASSICAL_SOLVER_CUH
# define CLASSICAL_SOLVER_CUH

# include "../eikonal.hpp"

class Classical : public Eikonal
{
private:

    int nit;

    float * K = nullptr;    
    float * d_K = nullptr;
    float * d_nK = nullptr;
    float * d_nT = nullptr;

    float * d_T = nullptr;
    float * d_S = nullptr;

    void set_volumes();
    void set_specifics();
    void initialization();

public:

    void set_forward_solver();
    void free_space();
};

__global__ void fdm_operators(float * S, float * T, float * K, float * nT, float h, int nxx, int nyy, int nzz);
__global__ void expanding_box(float * K, float * nK, int nxx, int nyy, int nzz);
__global__ void update_volume(float * T, float * nT, float * K, float * nK, int N);

# endif