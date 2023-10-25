# ifndef LEAST_SQUARES_CUH
# define LEAST_SQUARES_CUH

# include "../tomography.hpp"

class Least_Squares : public Tomography
{
private:

    int tk_order;   
    float lambda;

    int M, N, NNZ;

    int nx_tomo;  
    int ny_tomo;  
    int nz_tomo; 

    float dx_tomo;
    float dy_tomo;
    float dz_tomo;

    std::vector<int> iG;
    std::vector<int> jG;
    std::vector<float> vG;

    size_t ray_path_estimated_samples;

    int * iA = nullptr;
    int * jA = nullptr;
    float * vA = nullptr;
    float * B = nullptr;
    float * x = nullptr; 

    void apply_regularization();
    void solve_linear_system_lscg();
    void slowness_variation_rescaling();

    void apply_inversion_technique();
    void set_specific_parameters();
    void gradient_preconditioning();

public:

    void optimization();
};

# endif