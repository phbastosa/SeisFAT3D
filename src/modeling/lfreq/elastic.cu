# include "elastic.cuh"

void Elastic::set_specifications()
{
    fmax = std::stof(catch_parameter("max_frequency", parameters));

    set_wavelet();
    set_properties();    
    set_conditions();    

    nThreads = 256;
    nBlocks = (int)((volsize + nThreads - 1) / nThreads);

    current_xrec = new int[max_spread]();
    current_yrec = new int[max_spread]();
    current_zrec = new int[max_spread]();

    cudaMalloc((void**)&(rIdx), max_spread*sizeof(int));
    cudaMalloc((void**)&(rIdy), max_spread*sizeof(int));
    cudaMalloc((void**)&(rIdz), max_spread*sizeof(int));
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

void Elastic::set_wavelet()
{
    float * aux_s = new float[nt]();
    float * signal = new float[nt]();

    float pi = 4.0f*atanf(1.0f);
    float t0 = 2.0f*sqrtf(pi) / fmax;
    float fc = fmax / (3.0f * sqrtf(pi));

    tlag = (int)(t0 / dt) + 1;

    for (int n = 0; n < nt; n++)
    {
        float td = n*dt - t0;

        float arg = pi*pi*pi*fc*fc*td*td;

        aux_s[n] = 1e5f*(1.0f - 2.0f*arg)*expf(-arg);
    }

    for (int n = 0; n < nt; n++)
    {
        float summation = 0;
        for (int i = 0; i < n; i++)
            summation += aux_s[i];    
        
        signal[n] = summation;
    }

    cudaMalloc((void**)&(wavelet), nt*sizeof(float));

    cudaMemcpy(wavelet, signal, nt*sizeof(float), cudaMemcpyHostToDevice);

    delete[] signal;
}

void Elastic::export_synthetic_data()
{
    std::string data_file = data_folder + modeling_type + "_nStations" + std::to_string(spread) + "_nSamples" + std::to_string(nt) + "_shot_" + std::to_string(geometry->sInd[srcId]+1) + ".bin";
    export_binary_float(data_file, synthetic_data, nt*spread);    
}
