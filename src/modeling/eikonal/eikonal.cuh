# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include "../modeling.hpp"

class Eikonal : public Modeling
{
private:

    int iDivUp(int a, int b);

    void set_model_parameters();
    void set_model_boundaries();
    void set_wavefields();
    void set_outputs();

    void get_travelTimes();
    void get_firstArrivals();

protected:

    int totalLevels;

    int nSweeps, meshDim;

    float source_travel_time;    
    float dx2i, dy2i, dz2i, dsum; 
    float dx2dy2, dz2dx2, dz2dy2;

    float * S = nullptr;
    float * T = nullptr;
    float * V = nullptr;

    float * d_T = nullptr;
    float * d_S = nullptr;

    int * d_sgnv = nullptr;
    int * d_sgnt = nullptr;

public:

    void initial_setup();
    void forward_solver();
    void build_outputs();
    void free_space();
};

__global__ void fast_sweeping_kernel(float * S, float * T, int * sgnt, int * sgnv, int sgni, int sgnj, int sgnk, 
                                     int level, int xOffset, int yOffset, int xSweepOffset, int ySweepOffset, int zSweepOffset, 
                                     int nxx, int nyy, int nzz, float dx, float dy, float dz, float dx2i, float dy2i, float dz2i, 
                                     float dz2dx2, float dz2dy2, float dx2dy2, float dsum);


# endif