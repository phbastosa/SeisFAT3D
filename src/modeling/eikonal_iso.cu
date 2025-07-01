# include "eikonal_iso.cuh"

void Eikonal_ISO::set_conditions()
{
    modeling_type = "eikonal_iso";
    modeling_name = "Modeling type: Eikonal isotropic solver";
}

void Eikonal_ISO::time_propagation()
{
    initialization();
    eikonal_solver();
}
