# ifndef EIKONAL_VTI_CUH
# define EIKONAL_VTI_CUH

# include "eikonal.cuh"

class Eikonal_VTI : public Eikonal
{
private:

    int n, v, aId;

    float * E = nullptr;
    float * D = nullptr;
    float * G = nullptr;

    float * p = nullptr;
    float * Gv = nullptr;
    float * qV = nullptr;
    float * Gij = nullptr;
    float * Cijkl = nullptr;
    float * S_vti = nullptr;

    void get_stiffness();
    void get_christoffel();
    void get_eigen_values(); 

    void set_properties();
    void set_conditions();

    int voigt_map(int i, int j); 

public:

    void forward_solver();
};

# endif