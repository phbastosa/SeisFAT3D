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

    float * signal = new float[nt]();

    float pi = 4.0f*atanf(1.0f);

    float t0 = 2.0f*sqrtf(pi)/fmax;
    float fc = fmax/(3.0f * sqrtf(pi));

    for (int n = 0; n < nt; n++)
    {
        float td = n*dt - t0;

        float arg = pi*pi*pi*fc*fc*td*td;

        signal[n] = 1e5f*(1.0f - 2.0f*arg)*expf(-arg);
    }
    
    cudaMalloc((void**)&(wavelet), nt*sizeof(float));

    cudaMemcpy(wavelet, signal, nt*sizeof(float), cudaMemcpyHostToDevice);    

    cudaMalloc((void**)&(P), volsize*sizeof(float));
    cudaMalloc((void**)&(Pold), volsize*sizeof(float));
    cudaMalloc((void**)&(Pnew), volsize*sizeof(float));
}

void Scalar::initialization()
{
    cudaMemset(P, 0.0f, volsize*sizeof(float));
    cudaMemset(Pold, 0.0f, volsize*sizeof(float));
    cudaMemset(Pnew, 0.0f, volsize*sizeof(float));

    snap_index = 0;
}

void Scalar::set_forward_solver()
{
    for (time_index = 0; time_index < nt; time_index++)
    {
        display_progress();

        compute_pressure<<<blocksPerGrid,threadsPerBlock>>>(P,Pold,Pnew,V,damp1D,damp2D,damp3D,wavelet,source_index,time_index,dx,dy,dz,dt,nxx,nyy,nzz,nabc);
        cudaDeviceSynchronize();

        update_pressure<<<blocksPerGrid,threadsPerBlock>>>(P,Pold,Pnew,volsize);
        cudaDeviceSynchronize();

        get_wavefield_output();
        get_seismogram();
    }   

    get_receiver_output();
}

void Scalar::free_space()
{
    cudaFree(P);
    cudaFree(Pold);
    cudaFree(Pnew);
}

__global__ void compute_pressure(float * P, float * Pold, float * Pnew, float * V, float * damp1D, float * damp2D, float * damp3D, float * wavelet, int sId, int tId, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz, int nabc)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         // y direction
    int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
    int i = (int) (index - j*nzz - k*nxx*nzz); // z direction
    
    if (index == 0) P[sId] += wavelet[tId] / (dx*dy*dz);

    if((i >= 4) && (i < nzz-4) && (j >= 4) && (j < nxx-4) && (k >= 4) && (k < nyy-4)) 
    {
        float d2P_dx2 = (- 9.0f*(P[i + (j-4)*nzz + k*nxx*nzz] + P[i + (j+4)*nzz + k*nxx*nzz])
                     +   128.0f*(P[i + (j-3)*nzz + k*nxx*nzz] + P[i + (j+3)*nzz + k*nxx*nzz])
                     -  1008.0f*(P[i + (j-2)*nzz + k*nxx*nzz] + P[i + (j+2)*nzz + k*nxx*nzz])
                     +  8064.0f*(P[i + (j-1)*nzz + k*nxx*nzz] + P[i + (j+1)*nzz + k*nxx*nzz])
                     - 14350.0f*(P[i + j*nzz + k*nxx*nzz]))/(5040.0f*powf(dx, 2.0f));

        float d2P_dy2 = (- 9.0f*(P[i + j*nzz + (k-4)*nxx*nzz] + P[i + j*nzz + (k+4)*nxx*nzz])
                     +   128.0f*(P[i + j*nzz + (k-3)*nxx*nzz] + P[i + j*nzz + (k+3)*nxx*nzz])
                     -  1008.0f*(P[i + j*nzz + (k-2)*nxx*nzz] + P[i + j*nzz + (k+2)*nxx*nzz])
                     +  8064.0f*(P[i + j*nzz + (k-1)*nxx*nzz] + P[i + j*nzz + (k+1)*nxx*nzz])
                     - 14350.0f*(P[i + j*nzz + k*nxx*nzz]))/(5040.0f*powf(dy,2.0f));

        float d2P_dz2 = (- 9.0f*(P[(i-4) + j*nzz + k*nxx*nzz] + P[(i+4) + j*nzz + k*nxx*nzz])
                     +   128.0f*(P[(i-3) + j*nzz + k*nxx*nzz] + P[(i+3) + j*nzz + k*nxx*nzz])
                     -  1008.0f*(P[(i-2) + j*nzz + k*nxx*nzz] + P[(i+2) + j*nzz + k*nxx*nzz])
                     +  8064.0f*(P[(i-1) + j*nzz + k*nxx*nzz] + P[(i+1) + j*nzz + k*nxx*nzz])
                     - 14350.0f*(P[i + j*nzz + k*nxx*nzz]))/(5040.0f*powf(dz,2.0f));
    
        Pnew[index] = dt*dt*V[index]*V[index] * (d2P_dx2 + d2P_dy2 + d2P_dz2) + 2.0f*P[index] - Pold[index]; 
    }

    float damper = get_boundary_damper(damp1D,damp2D,damp3D,i,j,k,nxx,nyy,nzz,nabc);
    
    if (index < nxx*nyy*nzz)
    {
        P[index] *= damper;
        Pold[index] *= damper;    
        Pnew[index] *= damper;
    }
}

__global__ void update_pressure(float * P, float * Pold, float * Pnew, int volsize)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (index < volsize)
    {
        Pold[index] = P[index];        
        P[index] = Pnew[index];
    }
}