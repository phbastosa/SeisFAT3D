# include "acoustic.cuh"

void Acoustic::set_models()
{
    std::string vp_file = catch_parameter("vp_model_file", file);
    std::string rho_file = catch_parameter("rho_model_file", file);

    float * v = new float[nPoints]();
    float * p = new float[nPoints]();

    float * k = new float[volsize]();
    float * b = new float[volsize]();

    import_binary_float(vp_file, v, nPoints);
    import_binary_float(rho_file, p, nPoints);

    expand_boundary(v, k);
    expand_boundary(p, b);

    for (int index = 0; index < volsize; index++)
    {
        k[index] = b[index]*k[index]*k[index];
        b[index] = 1.0f / b[index];
    }

    cudaMalloc((void**)&(K), volsize*sizeof(float));
    cudaMalloc((void**)&(B), volsize*sizeof(float));
    
    cudaMemcpy(K, k, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(B, b, volsize*sizeof(float), cudaMemcpyHostToDevice);

    delete[] v;    
    delete[] p;    
    delete[] k;    
    delete[] b;    
}

void Acoustic::set_volumes()
{
    type_name = std::string("acoustic");
    type_message = std::string("[4] - Acoustic isotropic media");

    define_staggered_wavelet();

    cudaMalloc((void**)&(P), volsize*sizeof(float));
    cudaMalloc((void**)&(Vx), volsize*sizeof(float));
    cudaMalloc((void**)&(Vy), volsize*sizeof(float));
    cudaMalloc((void**)&(Vz), volsize*sizeof(float));
}

void Acoustic::initialization()
{
    cudaMemset(P, 0.0f, volsize*sizeof(float));
    cudaMemset(Vx, 0.0f, volsize*sizeof(float));
    cudaMemset(Vy, 0.0f, volsize*sizeof(float));
    cudaMemset(Vz, 0.0f, volsize*sizeof(float));
}

void Acoustic::set_forward_solver()
{
    for (time_index = 0; time_index < nt; time_index++)
    {
        // display_progression();

        // apply_wavelet<<<1,1>>>();

        // fdm_8E2T_get_velocity<<<>>>();
        // fdm_8E2T_get_pressure<<<>>>();

        // get_snapshots<<<>>>();
        // get_seismogram<<<>>>();
    }   
}