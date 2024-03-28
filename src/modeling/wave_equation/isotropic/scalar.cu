# include "scalar.cuh"

void Scalar::set_models()
{
    std::string vp_file = catch_parameter("vp_model_file", file);

    float * vp = new float[nPoints]();
    float * v = new float[volsize]();

    import_binary_float(vp_file, vp, nPoints);

    expand_boundary(vp, v);

    cudaMalloc((void**)&(V), volsize*sizeof(float));
    
    cudaMemcpy(V, v, volsize*sizeof(float), cudaMemcpyHostToDevice);

    delete[] v;
    delete[] vp;    
}

void Scalar::set_volumes()
{
    type_name = std::string("scalar");
    type_message = std::string("[3] - Scalar isotropic media");

    define_common_wavelet();

    cudaMalloc((void**)&(P), volsize*sizeof(float));
    cudaMalloc((void**)&(Pold), volsize*sizeof(float));
    cudaMalloc((void**)&(Pnew), volsize*sizeof(float));
}

void Scalar::initialization()
{
    cudaMemset(P, 0.0f, volsize*sizeof(float));
    cudaMemset(Pold, 0.0f, volsize*sizeof(float));
    cudaMemset(Pnew, 0.0f, volsize*sizeof(float));
}

void Scalar::set_forward_solver()
{
    for (time_index = 0; time_index < nt; time_index++)
    {
        // display_progression();

        // apply_wavelet<<<1,1>>>();

        // fdm_8E2T_scalar<<<>>>();

        // update_wavefield<<<>>>();    

        // get_snapshots<<<>>>();
        // get_seismogram<<<>>>();
    }   
}