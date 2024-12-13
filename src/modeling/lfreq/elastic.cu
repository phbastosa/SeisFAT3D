# include "elastic.cuh"

void Elastic::set_specifications()
{
    fmax = std::stof(catch_parameter("max_frequency", parameters));

    set_wavelet();
    set_boundaries();
    set_properties();    
    set_conditions();    

    nThreads = 256;
    nBlocks = (int)((volsize + nThreads - 1) / nThreads);

    P = new float[volsize]();

    current_xrec = new int[max_spread]();
    current_yrec = new int[max_spread]();
    current_zrec = new int[max_spread]();

    synthetic_data = new float[nt*max_spread]();

    cudaMalloc((void**)&(d_P), volsize*sizeof(float));

    cudaMalloc((void**)&(d_Vx), volsize*sizeof(float));
    cudaMalloc((void**)&(d_Vy), volsize*sizeof(float));
    cudaMalloc((void**)&(d_Vz), volsize*sizeof(float));

    cudaMalloc((void**)&(d_Txx), volsize*sizeof(float));
    cudaMalloc((void**)&(d_Tyy), volsize*sizeof(float));
    cudaMalloc((void**)&(d_Tzz), volsize*sizeof(float));
    
	cudaMalloc((void**)&(d_Txz), volsize*sizeof(float));
	cudaMalloc((void**)&(d_Tyz), volsize*sizeof(float));
	cudaMalloc((void**)&(d_Txy), volsize*sizeof(float));

    cudaMalloc((void**)&(rIdx), max_spread*sizeof(int));
    cudaMalloc((void**)&(rIdy), max_spread*sizeof(int));
    cudaMalloc((void**)&(rIdz), max_spread*sizeof(int));

    cudaMalloc((void**)&(seismogram), nt*max_spread*sizeof(float));
}

void Elastic::set_wavelet()
{
    float * signal_aux1 = new float[nt]();
    float * signal_aux2 = new float[nt]();

    float pi = 4.0f*atanf(1.0f);
    float t0 = 2.0f*sqrtf(pi) / fmax;
    float fc = fmax / (3.0f * sqrtf(pi));

    tlag = (int)(t0 / dt) + 1;

    for (int n = 0; n < nt; n++)
    {
        float td = n*dt - t0;

        float arg = pi*pi*pi*fc*fc*td*td;

        signal_aux1[n] = 1e5f*(1.0f - 2.0f*arg)*expf(-arg);
    }

    for (int n = 0; n < nt; n++)
    {
        float summation = 0;
        for (int i = 0; i < n; i++)
            summation += signal_aux1[i];    
        
        signal_aux2[n] = summation;
    }

    cudaMalloc((void**)&(wavelet), nt*sizeof(float));

    cudaMemcpy(wavelet, signal_aux2, nt*sizeof(float), cudaMemcpyHostToDevice);

    delete[] signal_aux1;
    delete[] signal_aux2;
}

void Elastic::set_boundaries()
{
    nb = std::stoi(catch_parameter("boundary_samples", parameters));
    bd = std::stof(catch_parameter("boundary_damping", parameters));

    float * damp1D = new float[nb]();
    float * damp2D = new float[nb*nb]();
    float * damp3D = new float[nb*nb*nb]();

    for (int i = 0; i < nb; i++) 
    {
        damp1D[i] = expf(-powf(bd * (nb - i), 2.0f));
    }

    for(int i = 0; i < nb; i++) 
    {
        for (int j = 0; j < nb; j++)
        {   
            damp2D[j + i*nb] += damp1D[i];
            damp2D[i + j*nb] += damp1D[i];
        }
    }

    for (int i  = 0; i < nb; i++)
    {
        for(int j = 0; j < nb; j++)
        {
            for(int k = 0; k < nb; k++)
            {
                damp3D[i + j*nb + k*nb*nb] += damp2D[i + j*nb];
                damp3D[i + j*nb + k*nb*nb] += damp2D[j + k*nb];
                damp3D[i + j*nb + k*nb*nb] += damp2D[i + k*nb];
            }
        }
    }    

    for (int index = 0; index < nb*nb; index++)
        damp2D[index] -= 1.0f;

    for (int index = 0; index < nb*nb*nb; index++)
        damp3D[index] -= 5.0f;    

	cudaMalloc((void**)&(d1D), nb*sizeof(float));
	cudaMalloc((void**)&(d2D), nb*nb*sizeof(float));
	cudaMalloc((void**)&(d3D), nb*nb*nb*sizeof(float));

	cudaMemcpy(d1D, damp1D, nb*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d2D, damp2D, nb*nb*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d3D, damp3D, nb*nb*nb*sizeof(float), cudaMemcpyHostToDevice);

    delete[] damp1D;
    delete[] damp2D;
    delete[] damp3D;
}

void Elastic::set_properties()
{
    std::string vp_file = catch_parameter("vp_model_file", parameters);
    std::string vs_file = catch_parameter("vs_model_file", parameters);
    std::string ro_file = catch_parameter("ro_model_file", parameters);

    float * vp = new float[nPoints]();
    float * vs = new float[nPoints]();
    float * ro = new float[nPoints]();

    Vp = new float[volsize]();
    Vs = new float[volsize]();
    Ro = new float[volsize]();

    import_binary_float(vp_file, vp, nPoints);
    import_binary_float(vs_file, vs, nPoints);
    import_binary_float(ro_file, ro, nPoints);

    expand_boundary(vp, Vp);
    expand_boundary(vs, Vs);
    expand_boundary(ro, Ro);

    delete[] vp;
    delete[] vs;
    delete[] ro;
}

void Elastic::initialization()
{
    cudaMemset(d_P, 0.0f, volsize*sizeof(float));
    
	cudaMemset(d_Vx, 0.0f, volsize*sizeof(float));
    cudaMemset(d_Vy, 0.0f, volsize*sizeof(float));
    cudaMemset(d_Vz, 0.0f, volsize*sizeof(float));
    
	cudaMemset(d_Txx, 0.0f, volsize*sizeof(float));
	cudaMemset(d_Tyy, 0.0f, volsize*sizeof(float));
    cudaMemset(d_Tzz, 0.0f, volsize*sizeof(float));
    
	cudaMemset(d_Txz, 0.0f, volsize*sizeof(float));
	cudaMemset(d_Tyz, 0.0f, volsize*sizeof(float));
	cudaMemset(d_Txy, 0.0f, volsize*sizeof(float));

    sIdx = (int)(geometry->xsrc[geometry->sInd[srcId]] / dx) + nb;
    sIdy = (int)(geometry->ysrc[geometry->sInd[srcId]] / dy) + nb;
    sIdz = (int)(geometry->zsrc[geometry->sInd[srcId]] / dz) + nb;

    int spread = 0;

    for (recId = geometry->iRec[srcId]; recId < geometry->fRec[srcId]; recId++)
    {
        current_xrec[spread] = (int)(geometry->xrec[recId] / dx) + nb;
        current_yrec[spread] = (int)(geometry->yrec[recId] / dy) + nb;
        current_zrec[spread] = (int)(geometry->zrec[recId] / dz) + nb;
    
        ++spread;
    }

    sBlocks = (int)((geometry->spread[srcId] + nThreads - 1) / nThreads); 

    cudaMemcpy(rIdx, current_xrec, geometry->spread[srcId]*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(rIdy, current_yrec, geometry->spread[srcId]*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(rIdz, current_zrec, geometry->spread[srcId]*sizeof(int), cudaMemcpyHostToDevice);

    cudaMemset(seismogram, 0.0f, nt*geometry->spread[srcId]*sizeof(float));
}

void Elastic::export_synthetic_data()
{
    std::string data_file = data_folder + modeling_type + "_nStations" + std::to_string(geometry->spread[srcId]) + "_nSamples" + std::to_string(nt) + "_shot_" + std::to_string(geometry->sInd[srcId]+1) + ".bin";
    export_binary_float(data_file, synthetic_data, nt*geometry->spread[srcId]);    
}
