# include "elastic_vti.cuh"

void Elastic_VTI::set_conditions()
{
    modeling_type = "elastic_vti";
    modeling_name = "Modeling type: Elastic vertically transverse isotropic time propagation";

    std::string e_file = catch_parameter("epsilon_model_file", parameters);
    std::string d_file = catch_parameter("delta_model_file", parameters);
    std::string g_file = catch_parameter("gamma_model_file", parameters);

    float * e = new float[nPoints]();
    float * d = new float[nPoints]();
    float * g = new float[nPoints]();

    E = new float[volsize]();
    D = new float[volsize]();
    G = new float[volsize]();

    import_binary_float(e_file, e, nPoints);
    import_binary_float(d_file, d, nPoints);
    import_binary_float(g_file, g, nPoints);

    expand_boundary(e, E);
    expand_boundary(d, D);
    expand_boundary(g, G);

    delete[] e;
    delete[] d;
    delete[] g;
}

void Elastic_VTI::forward_solver()
{


}

