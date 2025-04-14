# include "eikonal_iso.cuh"

void Eikonal_ISO::set_properties()
{
    Vp = new float[nPoints]();

    std::string model_file = catch_parameter("vp_model_file", parameters);

    import_binary_float(model_file, Vp, nPoints);

    for (int index = 0; index < nPoints; index++)
        Vp[index] = 1.0f / Vp[index];

    S = new float[volsize]();

    expand_boundary(Vp, S);
}

void Eikonal_ISO::set_conditions()
{
    modeling_type = "eikonal_iso";
    modeling_name = "Modeling type: Eikonal isotropic solver";
}

void Eikonal_ISO::forward_solver()
{
    initialization();

    propagation();

    compute_seismogram();
}
