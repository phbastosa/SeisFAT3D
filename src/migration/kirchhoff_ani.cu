# include "kirchhoff_ani.cuh"

void Kirchhoff_ANI::set_modeling_type()
{
    modeling = new Eikonal_ANI();
}