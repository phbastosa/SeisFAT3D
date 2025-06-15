# include "eikonal_iso.cuh"

void Eikonal_ISO::set_properties()
{
    float * vp = new float[nPoints]();

    std::string model_file = catch_parameter("vp_model_file", parameters);

    import_binary_float(model_file, vp, nPoints);

    # pragma omp parallel for
    for (int index = 0; index < nPoints; index++)
        vp[index] = 1.0f / vp[index];

    S = new float[volsize]();

    expand_boundary(vp, S);

    delete[] vp;
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
