# include "acoustic_isotropic.cuh"

void Acoustic_Isotropic::set_modeling_message()
{
    std::cout<<"Running:\n";
    std::cout<<"[4] - Solution for wave equation in variable density acoustic isotropic media\n\n"; 

    std::cout<<"Modeling progress: " << floorf(100.0f * (float)(time_id+1) / (float)(nt)) <<" %\n\n";
}

void Acoustic_Isotropic::set_model_parameters()
{
    modeling_method = std::string("acoustic_isotropic");

    float * vp = new float[nPoints]();
    float * rho = new float[nPoints]();

    // import_binary_float(catch_parameter("vp_model_file", file), vp, nPoints);
    // import_binary_float(catch_parameter("rho_model_file", file), rho, nPoints);

    for (int index = 0; index < nPoints; index++) 
    {
        vp[index] = 1500.0f;
        rho[index] = 1000.0f;
    }   

    Vp = new float[volsize]();
    Rho = new float[volsize]();

    expand_boundary(vp, Vp);
    expand_boundary(rho, Rho);

    delete[] vp;
    delete[] rho;

    float * b = new float[volsize]();
    float * k = new float[volsize]();

    for (int index = 0; index < volsize; index++)
    {
        b[index] = 1.0f / Rho[index];
        k[index] = Vp[index]*Vp[index]*Rho[index];
    }
    
	cudaMalloc((void**)&(B), volsize*sizeof(float));
	cudaMalloc((void**)&(K), volsize*sizeof(float));

	cudaMemcpy(B, b, volsize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(K, k, volsize*sizeof(float), cudaMemcpyHostToDevice);

    delete[] b;
    delete[] k;
}

void Acoustic_Isotropic::set_wavefields()
{
	cudaMalloc((void**)&(Vx), volsize*sizeof(float));
	cudaMalloc((void**)&(Vy), volsize*sizeof(float));
	cudaMalloc((void**)&(Vz), volsize*sizeof(float));

    cudaMalloc((void**)&(Pressure), volsize*sizeof(float));
}

void Acoustic_Isotropic::initial_setup()
{
    cudaMemset(Vx, 0.0f, volsize*sizeof(float));
    cudaMemset(Vy, 0.0f, volsize*sizeof(float));
    cudaMemset(Vz, 0.0f, volsize*sizeof(float));

    cudaMemset(Pressure, 0.0f, volsize*sizeof(float));
    cudaMemset(seismogram, 0.0f, nt*total_nodes*sizeof(float));

    int sidx = (int)(geometry->shots.x[shot_id] / dx) + nbxl;
    int sidy = (int)(geometry->shots.y[shot_id] / dy) + nbyl;
    int sidz = (int)(geometry->shots.z[shot_id] / dz) + nbzu;

    sId = sidz + sidx*nzz + sidy*nxx*nzz;

    isnap = 0;
}

void Acoustic_Isotropic::forward_solver()
{
    for (time_id = 0; time_id < nt; time_id++)
    {
        show_progress();
        
        // apply_wavelet

        // compute_velocity

        // compute_pressure
        
        get_snapshots();
        get_seismogram();
    }
}

void Acoustic_Isotropic::free_space()
{
    cudaFree(B);
    cudaFree(K);

    cudaFree(Vx);
    cudaFree(Vy);
    cudaFree(Vz);

    cudaFree(damp1D);
    cudaFree(damp2D);
    cudaFree(damp3D);
    
    cudaFree(wavelet);
    cudaFree(Pressure);    
    cudaFree(seismogram);
}
