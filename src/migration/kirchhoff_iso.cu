# include "kirchhoff_iso.cuh"

void Kirchhoff_ISO::set_modeling_type()
{
    modeling = new Eikonal_ISO();
}