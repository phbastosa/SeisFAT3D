# include "scalar.hpp"

void Scalar::set_wavelet()
{
    float * ricker = new float[nt]();

    for (int n = 0; n < nt; n++)
    {        
        float arg = pi*((n*dt - tlag)*fc*pi)*((n*dt - tlag)*fc*pi);
        
        ricker[n] = amp*(1 - 2*arg)*expf(-arg);    
    }

    if (import_wavelet) 
        import_binary_float(wavelet_file, ricker, nt);

	cudaMalloc((void**)&(wavelet), nt*sizeof(float));
	cudaMemcpy(wavelet, ricker, nt*sizeof(float), cudaMemcpyHostToDevice);

    delete[] ricker;
}