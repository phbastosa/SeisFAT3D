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

    snap_index = 0;
}

void Acoustic::set_forward_solver()
{
    for (time_index = 0; time_index < nt; time_index++)
    {
        display_progress();

        compute_velocity<<<blocksPerGrid,threadsPerBlock>>>(P,Vx,Vy,Vz,B,wavelet,source_index,time_index,dx,dy,dz,dt,nxx,nyy,nzz);
        cudaDeviceSynchronize();

        compute_pressure<<<blocksPerGrid,threadsPerBlock>>>(P,Vx,Vy,Vz,K,damp1D,damp2D,damp3D,dx,dy,dz,dt,nxx,nyy,nzz,nabc);
        cudaDeviceSynchronize();

        get_wavefield_output();
        get_seismogram();
    }   

    get_receiver_output();
}

void Acoustic::free_space()
{
    cudaFree(Vx);
    cudaFree(Vy);
    cudaFree(Vz);
    cudaFree(P);
}

__global__ void compute_velocity(float * P, float * Vx, float * Vy, float * Vz, float * B, float * wavelet, int sId, int tId, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         // y direction
    int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
    int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

    if (index == 0) P[sId] += wavelet[tId] / (dx*dy*dz);

    if((i > 3) && (i < nzz-3) && (j > 3) && (j < nxx-3) && (k > 3) && (k < nyy-3)) 
    {
        float dP_dx = (75.0f*(P[i + (j-4)*nzz + k*nxx*nzz] - P[i + (j+3)*nzz + k*nxx*nzz]) +
                     1029.0f*(P[i + (j+2)*nzz + k*nxx*nzz] - P[i + (j-3)*nzz + k*nxx*nzz]) +
                     8575.0f*(P[i + (j-2)*nzz + k*nxx*nzz] - P[i + (j+1)*nzz + k*nxx*nzz]) +
                   128625.0f*(P[i + j*nzz + k*nxx*nzz]     - P[i + (j-1)*nzz + k*nxx*nzz])) / (107520.0f*dx);

        float dP_dy = (75.0f*(P[i + j*nzz + (k-4)*nxx*nzz] - P[i + j*nzz + (k+3)*nxx*nzz]) +
                     1029.0f*(P[i + j*nzz + (k+2)*nxx*nzz] - P[i + j*nzz + (k-3)*nxx*nzz]) +
                     8575.0f*(P[i + j*nzz + (k-2)*nxx*nzz] - P[i + j*nzz + (k+1)*nxx*nzz]) +
                   128625.0f*(P[i + j*nzz + k*nxx*nzz]     - P[i + j*nzz + (k-1)*nxx*nzz])) / (107520.0f*dy);

        float dP_dz = (75.0f*(P[(i-4) + j*nzz + k*nxx*nzz] - P[(i+3) + j*nzz + k*nxx*nzz]) +
                     1029.0f*(P[(i+2) + j*nzz + k*nxx*nzz] - P[(i-3) + j*nzz + k*nxx*nzz]) +
                     8575.0f*(P[(i-2) + j*nzz + k*nxx*nzz] - P[(i+1) + j*nzz + k*nxx*nzz]) +
                   128625.0f*(P[i + j*nzz + k*nxx*nzz]     - P[(i-1) + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

        float Bx = 0.5f*(B[i + (j+1)*nzz + k*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);
        float By = 0.5f*(B[i + j*nzz + (k+1)*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);
        float Bz = 0.5f*(B[(i+1) + j*nzz + k*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);

        Vx[index] += -dt*Bx*dP_dx; 
        Vy[index] += -dt*By*dP_dy; 
        Vz[index] += -dt*Bz*dP_dz; 
    }
}

__global__ void compute_pressure(float * P, float * Vx, float * Vy, float * Vz, float * K, float * damp1D, float * damp2D, float * damp3D, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz, int nabc)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         // y direction
    int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
    int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

    if((i >= 3) && (i < nzz-4) && (j >= 3) && (j < nxx-4) && (k >= 3) && (k < nyy-4)) 
    {
        float dVx_dx = (75.0f*(Vx[i + (j-3)*nzz + k*nxx*nzz] - Vx[i + (j+4)*nzz + k*nxx*nzz]) +
                      1029.0f*(Vx[i + (j+3)*nzz + k*nxx*nzz] - Vx[i + (j-2)*nzz + k*nxx*nzz]) +
                      8575.0f*(Vx[i + (j-1)*nzz + k*nxx*nzz] - Vx[i + (j+2)*nzz + k*nxx*nzz]) +
                    128625.0f*(Vx[i + (j+1)*nzz + k*nxx*nzz] - Vx[i + j*nzz + k*nxx*nzz])) / (107520.0f*dx);

        float dVy_dy = (75.0f*(Vy[i + j*nzz + (k-3)*nxx*nzz] - Vy[i + j*nzz + (k+4)*nxx*nzz]) +
                      1029.0f*(Vy[i + j*nzz + (k+3)*nxx*nzz] - Vy[i + j*nzz + (k-2)*nxx*nzz]) +
                      8575.0f*(Vy[i + j*nzz + (k-1)*nxx*nzz] - Vy[i + j*nzz + (k+2)*nxx*nzz]) +
                    128625.0f*(Vy[i + j*nzz + (k+1)*nxx*nzz] - Vy[i + j*nzz + k*nxx*nzz])) / (107520.0f*dy);

        float dVz_dz = (75.0f*(Vz[(i-3) + j*nzz + k*nxx*nzz] - Vz[(i+4) + j*nzz + k*nxx*nzz]) +
                      1029.0f*(Vz[(i+3) + j*nzz + k*nxx*nzz] - Vz[(i-2) + j*nzz + k*nxx*nzz]) +
                      8575.0f*(Vz[(i-1) + j*nzz + k*nxx*nzz] - Vz[(i+2) + j*nzz + k*nxx*nzz]) +
                    128625.0f*(Vz[(i+1) + j*nzz + k*nxx*nzz] - Vz[i + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

        P[index] += -dt*K[index]*(dVx_dx + dVy_dy + dVz_dz); 
    }

    float damper = 1.0f;

    // 1D damping
    if((i < nabc) && (j >= nabc) && (j < nxx-nabc) && (k >= nabc) && (k < nyy-nabc)) 
    {
        damper = damp1D[i];
    }         
    else if((i >= nzz-nabc) && (i < nzz) && (j >= nabc) && (j < nxx-nabc) && (k >= nabc) && (k < nyy-nabc)) 
    {
        damper = damp1D[nabc-(i-(nzz-nabc))-1];
    }         
    else if((i >= nabc) && (i < nzz-nabc) && (j >= 0) && (j < nabc) && (k >= nabc) && (k < nyy-nabc)) 
    {
        damper = damp1D[j];
    }
    else if((i >= nabc) && (i < nzz-nabc) && (j >= nxx-nabc) && (j < nxx) && (k >= nabc) && (k < nyy-nabc)) 
    {
        damper = damp1D[nabc-(j-(nxx-nabc))-1];
    }
    else if((i >= nabc) && (i < nzz-nabc) && (j >= nabc) && (j < nxx-nabc) && (k >= 0) && (k < nabc)) 
    {
        damper = damp1D[k];
    }
    else if((i >= nabc) && (i < nzz-nabc) && (j >= nabc) && (j < nxx-nabc) && (k >= nyy-nabc) && (k < nyy)) 
    {
        damper = damp1D[nabc-(k-(nyy-nabc))-1];
    }

    // 2D damping 
    else if((i >= nabc) && (i < nzz-nabc) && (j >= 0) && (j < nabc) && (k >= 0) && (k < nabc))
    {
        damper = damp2D[j + k*nabc];
    }
    else if((i >= nabc) && (i < nzz-nabc) && (j >= nxx-nabc) && (j < nxx) && (k >= 0) && (k < nabc))
    {
        damper = damp2D[nabc-(j-(nxx-nabc))-1 + k*nabc];
    }
    else if((i >= nabc) && (i < nzz-nabc) && (j >= 0) && (j < nabc) && (k >= nyy-nabc) && (k < nyy))
    {
        damper = damp2D[j + (nabc-(k-(nyy-nabc))-1)*nabc];
    }
    else if((i >= nabc) && (i < nzz-nabc) && (j >= nxx-nabc) && (j < nxx) && (k >= nyy-nabc) && (k < nyy))
    {
        damper = damp2D[nabc-(j-(nxx-nabc))-1 + (nabc-(k-(nyy-nabc))-1)*nabc];
    }

    else if((i >= 0) && (i < nabc) && (j >= nabc) && (j < nxx-nabc) && (k >= 0) && (k < nabc))
    {
        damper = damp2D[i + k*nabc];
    }
    else if((i >= nzz-nabc) && (i < nzz) && (j >= nabc) && (j < nxx-nabc) && (k >= 0) && (k < nabc))
    {
        damper = damp2D[nabc-(i-(nzz-nabc))-1 + k*nabc];
    }
    else if((i >= 0) && (i < nabc) && (j >= nabc) && (j < nxx-nabc) && (k >= nyy-nabc) && (k < nyy))
    {
        damper = damp2D[i + (nabc-(k-(nyy-nabc))-1)*nabc];
    }
    else if((i >= nzz-nabc) && (i < nzz) && (j >= nabc) && (j < nxx-nabc) && (k >= nyy-nabc) && (k < nyy))
    {
        damper = damp2D[nabc-(i-(nzz-nabc))-1 + (nabc-(k-(nyy-nabc))-1)*nabc];
    }

    else if((i >= 0) && (i < nabc) && (j >= 0) && (j < nabc) && (k >= nabc) && (k < nyy-nabc))
    {
        damper = damp2D[i + j*nabc];
    }
    else if((i >= nzz-nabc) && (i < nzz) && (j >= 0) && (j < nabc) && (k >= nabc) && (k < nyy-nabc))
    {
        damper = damp2D[nabc-(i-(nzz-nabc))-1 + j*nabc];
    }
    else if((i >= 0) && (i < nabc) && (j >= nxx-nabc) && (j < nxx) && (k >= nabc) && (k < nyy-nabc))
    {
        damper = damp2D[i + (nabc-(j-(nxx-nabc))-1)*nabc];
    }
    else if((i >= nzz-nabc) && (i < nzz) && (j >= nxx-nabc) && (j < nxx) && (k >= nabc) && (k < nyy-nabc))
    {
        damper = damp2D[nabc-(i-(nzz-nabc))-1 + (nabc-(j-(nxx-nabc))-1)*nabc];
    }

    // 3D damping
    else if((i >= 0) && (i < nabc) && (j >= 0) && (j < nabc) && (k >= 0) && (k < nabc))
    {
        damper = damp3D[i + j*nabc + k*nabc*nabc];
    }
    else if((i >= nzz-nabc) && (i < nzz) && (j >= 0) && (j < nabc) && (k >= 0) && (k < nabc))
    {
        damper = damp3D[nabc-(i-(nzz-nabc))-1 + j*nabc + k*nabc*nabc];
    }
    else if((i >= 0) && (i < nabc) && (j >= nxx-nabc) && (j < nxx) && (k >= 0) && (k < nabc))
    {
        damper = damp3D[i + (nabc-(j-(nxx-nabc))-1)*nabc + k*nabc*nabc];
    }
    else if((i >= 0) && (i < nabc) && (j >= 0) && (j < nabc) && (k >= nyy-nabc) && (k < nyy))
    {
        damper = damp3D[i + j*nabc + (nabc-(k-(nyy-nabc))-1)*nabc*nabc];
    }
    else if((i >= nzz-nabc) && (i < nzz) && (j >= nxx-nabc) && (j < nxx) && (k >= 0) && (k < nabc))
    {
        damper = damp3D[nabc-(i-(nzz-nabc))-1 + (nabc-(j-(nxx-nabc))-1)*nabc + k*nabc*nabc];
    }
    else if((i >= nzz-nabc) && (i < nzz) && (j >= 0) && (j < nabc) && (k >= nyy-nabc) && (k < nyy))
    {
        damper = damp3D[nabc-(i-(nzz-nabc))-1 + j*nabc + (nabc-(k-(nyy-nabc))-1)*nabc*nabc];
    }
    else if((i >= 0) && (i < nabc) && (j >= nxx-nabc) && (j < nxx) && (k >= nyy-nabc) && (k < nyy))
    {
        damper = damp3D[i + (nabc-(j-(nxx-nabc))-1)*nabc + (nabc-(k-(nyy-nabc))-1)*nabc*nabc];
    }
    else if((i >= nzz-nabc) && (i < nzz) && (j >= nxx-nabc) && (j < nxx) && (k >= nyy-nabc) && (k < nyy))
    {
        damper = damp3D[nabc-(i-(nzz-nabc))-1 + (nabc-(j-(nxx-nabc))-1)*nabc + (nabc-(k-(nyy-nabc))-1)*nabc*nabc];
    }

    if (index < nxx*nyy*nzz)
    {
        Vx[index] *= damper;
        Vy[index] *= damper;
        Vz[index] *= damper;
        
        P[index] *= damper;
    }   
}
