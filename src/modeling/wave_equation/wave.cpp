# include "wave.hpp"

void Wave::set_specifics()
{
    nt = std::stoi(catch_parameter("time_samples", file));
    dt = std::stof(catch_parameter("time_spacing", file));

    fmax = std::stof(catch_parameter("max_frequency", file));

    nabc = std::stoi(catch_parameter("boundary_samples", file));
    pabc = std::stof(catch_parameter("boundary_damper", file));

    total_snaps = std::stoi(catch_parameter("total_snapshots", file));

    nbxl = nabc; nbxr = nabc;
    nbyl = nabc; nbyr = nabc;
    nbzu = nabc; nbzd = nabc;    

    define_cerjan_dampers();
    define_seismogram();
    define_snapshots();
}

void Wave::define_cerjan_dampers()
{
    float * d1D = new float[nabc]();
    float * d2D = new float[nabc*nabc]();
    float * d3D = new float[nabc*nabc*nabc]();

    float factor = std::stof(catch_parameter("boundary_damper", file));

    for (int i = 0; i < nabc; i++) 
    {
        d1D[i] = expf(-powf(factor * (nabc - i), 2.0f));
    }

    for(int i = 0; i < nabc; i++) 
    {
        for (int j = 0; j < nabc; j++)
        {   
            d2D[j + i*nabc] += d1D[i]; // up to bottom
            d2D[i + j*nabc] += d1D[i]; // left to right
        }
    }

    for (int i  = 0; i < nabc; i++)
    {
        for(int j = 0; j < nabc; j++)
        {
            for(int k = 0; k < nabc; k++)
            {
                d3D[i + j*nabc + k*nabc*nabc] += d2D[i + j*nabc]; // XY plane
                d3D[i + j*nabc + k*nabc*nabc] += d2D[j + k*nabc]; // ZX plane
                d3D[i + j*nabc + k*nabc*nabc] += d2D[i + k*nabc]; // ZY plane
            }
        }
    }    

    for (int index = 0; index < nabc*nabc; index++)
        d2D[index] -= 1.0f;

    for (int index = 0; index < nabc*nabc*nabc; index++)
        d3D[index] -= 5.0f;    

	cudaMalloc((void**)&(damp1D), nabc*sizeof(float));
	cudaMalloc((void**)&(damp2D), nabc*nabc*sizeof(float));
	cudaMalloc((void**)&(damp3D), nabc*nabc*nabc*sizeof(float));

	cudaMemcpy(damp1D, d1D, nabc*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(damp2D, d2D, nabc*nabc*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(damp3D, d3D, nabc*nabc*nabc*sizeof(float), cudaMemcpyHostToDevice);

    delete[] d1D;
    delete[] d2D;
    delete[] d3D;    
}

void Wave::define_seismogram()
{
    receiver_output_samples = nt*total_nodes;

    float * receiver_output = new float[receiver_output_samples]();

    cudaMalloc((void**)&(seismogram), receiver_output_samples*sizeof(float));
}

void Wave::define_snapshots()
{
    wavefield_output_samples = nPoints*total_snaps;

    float * wavefield_output = new float[wavefield_output_samples]();

    cudaMalloc((void**)&(snapshots), wavefield_output_samples*sizeof(float));
}

void Wave::define_common_wavelet()
{
    float * signal = new float[nt]();

    float pi = 4.0f*atanf(1.0f);

    float t0 = 2.0f*sqrtf(pi)/fmax;
    float fc = fmax/(3.0f * sqrtf(pi));

    for (int n = 0; n < nt; n++)
    {
        float td = n*dt - t0;

        float arg = pi*pi*fmax*fmax*td*td;

        signal[n] = (1.0f - 2.0f*pi*arg)*expf(-pi*arg);
    }
    
    cudaMalloc((void**)&(wavelet), nt*sizeof(float));

    cudaMemcpy(wavelet, signal, nt*sizeof(float), cudaMemcpyHostToDevice);

    delete[] signal;
}

void Wave::define_staggered_wavelet()
{
    float * signal = new float[nt]();

    float pi = 4.0f*atanf(1.0f);

    float t0 = 2.0f*sqrtf(pi)/fmax;
    float fc = fmax/(3.0f * sqrtf(pi));

    float summation = 0;

    for (int n = 0; n < nt; n++)
    {
        float td = n*dt - t0;

        float arg = pi*pi*fmax*fmax*td*td;

        summation += (1.0f - 2.0f*pi*arg)*expf(-pi*arg);

        signal[n] = summation;
    }
    
    cudaMalloc((void**)&(wavelet), nt*sizeof(float));

    cudaMemcpy(wavelet, signal, nt*sizeof(float), cudaMemcpyHostToDevice);

    delete[] signal;
}