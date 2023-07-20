# include "acoustic.hpp"

void Acoustic::set_wavelet()
{
    float * ricker = new float[nt]();
    float * aux = new float[nt]();

    float arg, sum;
    for (int n = 0; n < nt; n++)
    {        
        arg = pi*((n*dt - tlag)*fc*pi)*((n*dt - tlag)*fc*pi);
        
        ricker[n] = (1 - 2*arg) * expf(-arg);    

        sum = 0.0f;
        for (int k = 0; k < n+1; k++)     
            sum += ricker[k];
    
        aux[n] = amp * sum;
    }

	cudaMalloc((void**)&(wavelet), nt*sizeof(float));
	cudaMemcpy(wavelet, aux, nt*sizeof(float), cudaMemcpyHostToDevice);

    delete[] aux;
    delete[] ricker;
}

