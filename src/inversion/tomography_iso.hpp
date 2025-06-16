# ifndef TOMOGRAPHY_ISO_HPP
# define TOMOGRAPHY_ISO_HPP

# include "inversion.hpp"

class Tomography_ISO : public Inversion
{
private:

    int tk_order;
    float tk_param;

    size_t ray_path_max_samples;

    std::vector<int> iG;
    std::vector<int> jG;
    std::vector<float> vG;

    int M, N, NNZ;

    int * iA = nullptr;
    int * jA = nullptr;
    float * vA = nullptr;
    float * B = nullptr;
    float * x = nullptr; 

    void set_specifications();

    void apply_regularization();
    void solve_linear_system_lscg();
    void apply_inversion_technique();

public:

    void optimization();
};

# endif