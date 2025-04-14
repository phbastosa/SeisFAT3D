# ifndef EIKONAL_ANI_CUH
# define EIKONAL_ANI_CUH

# include "eikonal.cuh"

class Eikonal_ANI : public Eikonal
{
private:

    int n, v, aId;

    uintc * C11 = nullptr; float maxC11; float minC11;
    uintc * C12 = nullptr; float maxC12; float minC12;
    uintc * C13 = nullptr; float maxC13; float minC13;
    uintc * C14 = nullptr; float maxC14; float minC14;
    uintc * C15 = nullptr; float maxC15; float minC15;
    uintc * C16 = nullptr; float maxC16; float minC16;

    uintc * C22 = nullptr; float maxC22; float minC22;
    uintc * C23 = nullptr; float maxC23; float minC23;
    uintc * C24 = nullptr; float maxC24; float minC24;
    uintc * C25 = nullptr; float maxC25; float minC25;
    uintc * C26 = nullptr; float maxC26; float minC26;

    uintc * C33 = nullptr; float maxC33; float minC33;
    uintc * C34 = nullptr; float maxC34; float minC34;
    uintc * C35 = nullptr; float maxC35; float minC35;
    uintc * C36 = nullptr; float maxC36; float minC36;

    uintc * C44 = nullptr; float maxC44; float minC44;
    uintc * C45 = nullptr; float maxC45; float minC45;
    uintc * C46 = nullptr; float maxC46; float minC46;

    uintc * C55 = nullptr; float maxC55; float minC55;
    uintc * C56 = nullptr; float maxC56; float minC56;

    uintc * C66 = nullptr; float maxC66; float minC66;

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