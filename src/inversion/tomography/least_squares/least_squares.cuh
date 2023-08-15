# ifndef LEAST_SQUARES_CUH
# define LEAST_SQUARES_CUH

# include "../tomography.hpp"

class Least_Squares : public Tomography
{
private:

    int nx_tomo;  
    int ny_tomo;  
    int nz_tomo; 

    float dx_tomo;
    float dy_tomo;
    float dz_tomo;

    std::vector<int> iG;
    std::vector<int> jG;
    std::vector<float> vG;

    int M, N, NNZ;

    int tk_order;   
    float lambda;

    int * iA = nullptr;
    int * jA = nullptr;
    float * vA = nullptr;
    float * B = nullptr;
    float * x = nullptr;             // A x = B

    void gradient_ray_tracing();

public:

    void optimization();
    void set_parameters();
    void forward_modeling();
};

# endif