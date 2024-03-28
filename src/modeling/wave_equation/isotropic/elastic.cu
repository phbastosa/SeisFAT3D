# include "elastic.cuh"

void Elastic::set_models()
{
    std::string vp_file = catch_parameter("vp_model_file", file);
    std::string vs_file = catch_parameter("vs_model_file", file);
    std::string rho_file = catch_parameter("rho_model_file", file);

    float * vp = new float[nPoints]();
    float * vs = new float[nPoints]();
    float * p = new float[nPoints]();

    import_binary_float(vp_file, vp, nPoints);
    import_binary_float(vs_file, vs, nPoints);
    import_binary_float(rho_file, p, nPoints);

    float * b = new float[volsize]();
    float * l = new float[volsize]();
    float * m = new float[volsize]();

    expand_boundary(vp, l);
    expand_boundary(vs, m);
    expand_boundary(p, b);

    for (int index = 0; index < volsize; index++)
    {
        m[index] = b[index]*m[index]*m[index];
        l[index] = b[index]*l[index]*l[index] - 2.0f*m[index];
        b[index] = 1.0f / b[index];
    }

    cudaMalloc((void**)&(L), volsize*sizeof(float));
    cudaMalloc((void**)&(M), volsize*sizeof(float));
    cudaMalloc((void**)&(B), volsize*sizeof(float));
    
    cudaMemcpy(L, l, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(M, m, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(B, b, volsize*sizeof(float), cudaMemcpyHostToDevice);

    delete[] vp;    
    delete[] vs;    
    delete[] p;

    delete[] l;    
    delete[] m;    
    delete[] b;    
}

void Elastic::set_volumes()
{
    type_name = std::string("elastic");
    type_message = std::string("[5] - Elastic isotropic media");

    define_staggered_wavelet();

    cudaMalloc((void**)&(P), volsize*sizeof(float));
    cudaMalloc((void**)&(Vx), volsize*sizeof(float));
    cudaMalloc((void**)&(Vy), volsize*sizeof(float));
    cudaMalloc((void**)&(Vz), volsize*sizeof(float));
    cudaMalloc((void**)&(Txx), volsize*sizeof(float));
    cudaMalloc((void**)&(Tyy), volsize*sizeof(float));
    cudaMalloc((void**)&(Tzz), volsize*sizeof(float));
    cudaMalloc((void**)&(Txy), volsize*sizeof(float));
    cudaMalloc((void**)&(Txz), volsize*sizeof(float));
    cudaMalloc((void**)&(Tyz), volsize*sizeof(float));
}

void Elastic::initialization()
{
    cudaMemset(P, 0.0f, volsize*sizeof(float));
    cudaMemset(Vx, 0.0f, volsize*sizeof(float));
    cudaMemset(Vy, 0.0f, volsize*sizeof(float));
    cudaMemset(Vz, 0.0f, volsize*sizeof(float));
    cudaMemset(Txx, 0.0f, volsize*sizeof(float));
    cudaMemset(Tyy, 0.0f, volsize*sizeof(float));
    cudaMemset(Tzz, 0.0f, volsize*sizeof(float));
    cudaMemset(Txy, 0.0f, volsize*sizeof(float));
    cudaMemset(Txz, 0.0f, volsize*sizeof(float));
    cudaMemset(Tyz, 0.0f, volsize*sizeof(float));
}

void Elastic::set_forward_solver()
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