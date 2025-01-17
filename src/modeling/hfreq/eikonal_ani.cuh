# ifndef EIKONAL_ANI_CUH
# define EIKONAL_ANI_CUH

# include "eikonal.cuh"

class Eikonal_ANI : public Eikonal
{
private:

    int n, v, aId;

    float * C11 = nullptr;
    float * C12 = nullptr;
    float * C13 = nullptr;
    float * C14 = nullptr;
    float * C15 = nullptr;
    float * C16 = nullptr;

    float * C22 = nullptr;
    float * C23 = nullptr;
    float * C24 = nullptr;
    float * C25 = nullptr;
    float * C26 = nullptr;

    float * C33 = nullptr;
    float * C34 = nullptr;
    float * C35 = nullptr;
    float * C36 = nullptr;

    float * C44 = nullptr;
    float * C45 = nullptr;
    float * C46 = nullptr;

    float * C55 = nullptr;
    float * C56 = nullptr;

    float * C66 = nullptr;

    float * p = nullptr;
    float * C = nullptr;
    float * G = nullptr;
    float * Gv = nullptr;
    float * qS = nullptr;

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