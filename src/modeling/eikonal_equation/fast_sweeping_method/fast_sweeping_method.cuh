# ifndef FAST_SWEEPING_METHOD_HPP
# define FAST_SWEEPING_METHOD_HPP

# include "../eikonal.hpp"

class Fast_Sweeping_Method : public Eikonal
{
private:
        

protected:

    int totalLevels;
    int nSweeps, meshDim;

    float dx2i, dy2i, dz2i, dsum; 
    float dx2dy2, dz2dx2, dz2dy2;

    float * d_T = nullptr;
    float * d_S = nullptr;

    int * d_sgnv = nullptr;
    int * d_sgnt = nullptr;

    int iDivUp(int a, int b);

    void set_modeling_message();
    void set_model_boundaries();
    void set_preconditioners();

public:

    void initial_setup();
    void forward_solver();
    void free_space();
};

__global__ void fast_sweeping_kernel(float * S, float * T, int * sgnt, int * sgnv, int sgni, int sgnj, int sgnk, 
                                     int level, int xOffset, int yOffset, int xSweepOffset, int ySweepOffset, int zSweepOffset, 
                                     int nxx, int nyy, int nzz, float dx, float dy, float dz, float dx2i, float dy2i, float dz2i, 
                                     float dz2dx2, float dz2dy2, float dx2dy2, float dsum);


# endif