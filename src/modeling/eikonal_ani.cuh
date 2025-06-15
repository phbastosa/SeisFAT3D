# ifndef EIKONAL_ANI_CUH
# define EIKONAL_ANI_CUH

# include "eikonal.cuh"

class Eikonal_ANI : public Eikonal
{
private:

    float * qS = nullptr;

    uintc * d_C11 = nullptr; float maxC11; float minC11;
    uintc * d_C12 = nullptr; float maxC12; float minC12;
    uintc * d_C13 = nullptr; float maxC13; float minC13;
    uintc * d_C14 = nullptr; float maxC14; float minC14;
    uintc * d_C15 = nullptr; float maxC15; float minC15;
    uintc * d_C16 = nullptr; float maxC16; float minC16;

    uintc * d_C22 = nullptr; float maxC22; float minC22;
    uintc * d_C23 = nullptr; float maxC23; float minC23;
    uintc * d_C24 = nullptr; float maxC24; float minC24;
    uintc * d_C25 = nullptr; float maxC25; float minC25;
    uintc * d_C26 = nullptr; float maxC26; float minC26;

    uintc * d_C33 = nullptr; float maxC33; float minC33;
    uintc * d_C34 = nullptr; float maxC34; float minC34;
    uintc * d_C35 = nullptr; float maxC35; float minC35;
    uintc * d_C36 = nullptr; float maxC36; float minC36;

    uintc * d_C44 = nullptr; float maxC44; float minC44;
    uintc * d_C45 = nullptr; float maxC45; float minC45;
    uintc * d_C46 = nullptr; float maxC46; float minC46;

    uintc * d_C55 = nullptr; float maxC55; float minC55;
    uintc * d_C56 = nullptr; float maxC56; float minC56;

    uintc * d_C66 = nullptr; float maxC66; float minC66;

    void set_properties();
    void set_conditions();
 
public:

    void forward_solver();
};

__global__ void get_quasi_slowness(float * T, float * S, float dx, float dy, float dz, int sIdx, int sIdy, int sIdz, int nxx, int nyy, int nzz, int nb, 
                                   uintc * C11, uintc * C12, uintc * C13, uintc * C14, uintc * C15, uintc * C16, uintc * C22, uintc * C23, uintc * C24, uintc * C25, 
                                   uintc * C26, uintc * C33, uintc * C34, uintc * C35, uintc * C36, uintc * C44, uintc * C45, uintc * C46, uintc * C55, uintc * C56, 
                                   uintc * C66, float minC11, float maxC11, float minC12, float maxC12, float minC13, float maxC13, float minC14, float maxC14, 
                                   float minC15, float maxC15, float minC16, float maxC16, float minC22, float maxC22, float minC23, float maxC23, float minC24, 
                                   float maxC24, float minC25, float maxC25, float minC26, float maxC26, float minC33, float maxC33, float minC34, float maxC34, 
                                   float minC35, float maxC35, float minC36, float maxC36, float minC44, float maxC44, float minC45, float maxC45, float minC46, 
                                   float maxC46, float minC55, float maxC55, float minC56, float maxC56, float minC66, float maxC66);

# endif