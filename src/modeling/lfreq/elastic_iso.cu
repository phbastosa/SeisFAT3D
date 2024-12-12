# include "elastic_iso.cuh"

void Elastic_Iso::set_properties()
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

void Elastic_Iso::set_conditions()
{
    modeling_type = "elastic_iso";
    modeling_name = "Modeling type: Elastic isotropic wave propagation";

    M = new float[volsize]();
    L = new float[volsize]();
    B = new float[volsize]();
    P = new float[volsize]();

    for (int index = 0; index < volsize; index++)
    {
        M[index] = Ro[index]*Vs[index]*Vs[index];
        L[index] = Ro[index]*Vp[index]*Vp[index] - 2.0f*M[index];
        B[index] = 1.0f / Ro[index];
    }

    synthetic_data = new float[nt*max_spread]();
    cudaMalloc((void**)&(seismogram), nt*max_spread*sizeof(float));

    cudaMalloc((void**)&(d_M), volsize*sizeof(float));
    cudaMalloc((void**)&(d_L), volsize*sizeof(float));
    cudaMalloc((void**)&(d_B), volsize*sizeof(float));
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

    cudaMemcpy(d_M, M, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_L, L, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, B, volsize*sizeof(float), cudaMemcpyHostToDevice);
}

void Elastic_Iso::initialization()
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

void Elastic_Iso::forward_solver()
{
    for (int tId = 0; tId < nt + tlag; tId++)
    {
        compute_pressure<<<nBlocks, nThreads>>>(d_Vx, d_Vy, d_Vz, d_Txx, d_Tyy, d_Tzz, d_Txz, d_Tyz, d_Txy, d_P, d_M, d_L, wavelet, sIdx, sIdy, sIdz, tId, nt, dx, dy, dz, dt, nxx, nyy, nzz);
        cudaDeviceSynchronize();

        compute_velocity<<<nBlocks, nThreads>>>(d_Vx, d_Vy, d_Vz, d_Txx, d_Tyy, d_Tzz, d_Txz, d_Tyz, d_Txy, d_B, d1D, d2D, d3D, dx, dy, dz, dt, nxx, nyy, nzz, nb);
        cudaDeviceSynchronize();

        compute_seismogram<<<sBlocks, nThreads>>>(d_P, rIdx, rIdy, rIdz, seismogram, geometry->spread[srcId], tId, tlag, nt, nxx, nzz);     
        cudaDeviceSynchronize();
    }

    cudaMemcpy(synthetic_data, seismogram, nt*geometry->spread[srcId]*sizeof(float), cudaMemcpyDeviceToHost);
}

__global__ void compute_pressure(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * P, float * M, float * L, float * wavelet, int sIdx, int sIdy, int sIdz, int tId, int nt, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         // y direction
    int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
    int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

    if ((index == 0) && (tId < nt))
    {
        Txx[sIdz + sIdx*nzz + sIdy*nxx*nzz] += wavelet[tId] / (dx*dy*dz);
        Tyy[sIdz + sIdx*nzz + sIdy*nxx*nzz] += wavelet[tId] / (dx*dy*dz);
        Tzz[sIdz + sIdx*nzz + sIdy*nxx*nzz] += wavelet[tId] / (dx*dy*dz);
    }

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

        Txx[index] += dt*((L[index] + 2*M[index])*dVx_dx + L[index]*(dVy_dy + dVz_dz));
        Tyy[index] += dt*((L[index] + 2*M[index])*dVy_dy + L[index]*(dVx_dx + dVz_dz));
        Tzz[index] += dt*((L[index] + 2*M[index])*dVz_dz + L[index]*(dVx_dx + dVy_dy));                    
    }

    if((i >= 0) && (i < nzz) && (j > 3) && (j < nxx-3) && (k > 3) && (k < nyy-3)) 
    {
        float dVx_dy = (75.0f*(Vx[i + j*nzz + (k-4)*nxx*nzz] - Vx[i + j*nzz + (k+3)*nxx*nzz]) +
                      1029.0f*(Vx[i + j*nzz + (k+2)*nxx*nzz] - Vx[i + j*nzz + (k-3)*nxx*nzz]) +
                      8575.0f*(Vx[i + j*nzz + (k-2)*nxx*nzz] - Vx[i + j*nzz + (k+1)*nxx*nzz]) +
                    128625.0f*(Vx[i + j*nzz + k*nxx*nzz]     - Vx[i + j*nzz + (k-1)*nxx*nzz])) / (107520.0f*dy);

        float dVy_dx = (75.0f*(Vy[i + (j-4)*nzz + k*nxx*nzz] - Vy[i + (j+3)*nzz + k*nxx*nzz]) +
                      1029.0f*(Vy[i + (j+2)*nzz + k*nxx*nzz] - Vy[i + (j-3)*nzz + k*nxx*nzz]) +
                      8575.0f*(Vy[i + (j-2)*nzz + k*nxx*nzz] - Vy[i + (j+1)*nzz + k*nxx*nzz]) +
                    128625.0f*(Vy[i + j*nzz + k*nxx*nzz]     - Vy[i + (j-1)*nzz + k*nxx*nzz])) / (107520.0f*dx);

        float Mxy = powf(0.25f*(1.0f/M[i + (j+1)*nzz + (k+1)*nxx*nzz] + 1.0f/M[i + (j+1)*nzz + k*nxx*nzz] + 
                                1.0f/M[i + j*nzz + (k+1)*nxx*nzz]     + 1.0f/M[i + j*nzz + k*nxx*nzz]), -1.0f);

        Txy[index] += dt*Mxy*(dVx_dy + dVy_dx);
    }

    if((i > 3) && (i < nzz-3) && (j > 3) && (j < nxx-3) && (k >= 0) && (k < nyy)) 
    {
        float dVx_dz = (75.0f*(Vx[(i-4) + j*nzz + k*nxx*nzz] - Vx[(i+3) + j*nzz + k*nxx*nzz]) +
                      1029.0f*(Vx[(i+2) + j*nzz + k*nxx*nzz] - Vx[(i-3) + j*nzz + k*nxx*nzz]) +
                      8575.0f*(Vx[(i-2) + j*nzz + k*nxx*nzz] - Vx[(i+1) + j*nzz + k*nxx*nzz]) +
                    128625.0f*(Vx[i + j*nzz + k*nxx*nzz]     - Vx[(i-1) + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

        float dVz_dx = (75.0f*(Vz[i + (j-4)*nzz + k*nxx*nzz] - Vz[i + (j+3)*nzz + k*nxx*nzz]) +
                      1029.0f*(Vz[i + (j+2)*nzz + k*nxx*nzz] - Vz[i + (j-3)*nzz + k*nxx*nzz]) +
                      8575.0f*(Vz[i + (j-2)*nzz + k*nxx*nzz] - Vz[i + (j+1)*nzz + k*nxx*nzz]) +
                    128625.0f*(Vz[i + j*nzz + k*nxx*nzz]     - Vz[i + (j-1)*nzz + k*nxx*nzz])) / (107520.0f*dx);

        float Mxz = powf(0.25f*(1.0f/M[(i+1) + (j+1)*nzz + k*nxx*nzz] + 1.0f/M[i + (j+1)*nzz + k*nxx*nzz] + 
                                1.0f/M[(i+1) + j*nzz + k*nxx*nzz]     + 1.0f/M[i + j*nzz + k*nxx*nzz]), -1.0f);

        Txz[index] += dt*Mxz*(dVx_dz + dVz_dx);
    }

    if((i > 3) && (i < nzz-3) && (j >= 0) && (j < nxx) && (k > 3) && (k < nyy-3)) 
    {
        float dVy_dz = (75.0f*(Vy[(i-4) + j*nzz + k*nxx*nzz] - Vy[(i+3) + j*nzz + k*nxx*nzz]) +
                      1029.0f*(Vy[(i+2) + j*nzz + k*nxx*nzz] - Vy[(i-3) + j*nzz + k*nxx*nzz]) +
                      8575.0f*(Vy[(i-2) + j*nzz + k*nxx*nzz] - Vy[(i+1) + j*nzz + k*nxx*nzz]) +
                    128625.0f*(Vy[i + j*nzz + k*nxx*nzz]     - Vy[(i-1) + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

        float dVz_dy = (75.0f*(Vz[i + j*nzz + (k-4)*nxx*nzz] - Vz[i + j*nzz + (k+3)*nxx*nzz]) +
                      1029.0f*(Vz[i + j*nzz + (k+2)*nxx*nzz] - Vz[i + j*nzz + (k-3)*nxx*nzz]) +
                      8575.0f*(Vz[i + j*nzz + (k-2)*nxx*nzz] - Vz[i + j*nzz + (k+1)*nxx*nzz]) +
                    128625.0f*(Vz[i + j*nzz + k*nxx*nzz]     - Vz[i + j*nzz + (k-1)*nxx*nzz])) / (107520.0f*dy);

        float Myz = powf(0.25f*(1.0f/M[(i+1) + j*nzz + (k+1)*nxx*nzz] + 1.0f/M[i + j*nzz + (k+1)*nxx*nzz] + 
                                1.0f/M[(i+1) + j*nzz + k*nxx*nzz] +     1.0f/M[i + j*nzz + k*nxx*nzz]), -1.0f);

        Tyz[index] += dt*Myz*(dVy_dz + dVz_dy);
    }

    if ((i > 3) && (i < nzz-4) && (j > 3) && (j < nxx-4) && (k > 3) && (k < nyy-4))
    {
        P[index] = (Txx[index] + Tyy[index] + Tzz[index]) / 3.0f;
    }
}

__global__ void compute_velocity(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * B,  float * damp1D, float * damp2D, float * damp3D, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz, int nb)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         
    int j = (int) (index - k*nxx*nzz) / nzz;   
    int i = (int) (index - j*nzz - k*nxx*nzz); 

    if((i >= 3) && (i < nzz-4) && (j > 3) && (j < nxx-3) && (k >= 3) && (k < nyy-4)) 
    {
        float dTxx_dx = (75.0f*(Txx[i + (j-4)*nzz + k*nxx*nzz] - Txx[i + (j+3)*nzz + k*nxx*nzz]) +
                       1029.0f*(Txx[i + (j+2)*nzz + k*nxx*nzz] - Txx[i + (j-3)*nzz + k*nxx*nzz]) +
                       8575.0f*(Txx[i + (j-2)*nzz + k*nxx*nzz] - Txx[i + (j+1)*nzz + k*nxx*nzz]) +
                     128625.0f*(Txx[i + j*nzz + k*nxx*nzz]     - Txx[i + (j-1)*nzz + k*nxx*nzz])) / (107520.0f*dx);

        float dTxy_dy = (75.0f*(Txy[i + j*nzz + (k-3)*nxx*nzz] - Txy[i + j*nzz + (k+4)*nxx*nzz]) +
                       1029.0f*(Txy[i + j*nzz + (k+3)*nxx*nzz] - Txy[i + j*nzz + (k-2)*nxx*nzz]) +
                       8575.0f*(Txy[i + j*nzz + (k-1)*nxx*nzz] - Txy[i + j*nzz + (k+2)*nxx*nzz]) +
                     128625.0f*(Txy[i + j*nzz + (k+1)*nxx*nzz] - Txy[i + j*nzz + k*nxx*nzz])) / (107520.0f*dy);

        float dTxz_dz = (75.0f*(Txz[(i-3) + j*nzz + k*nxx*nzz] - Txz[(i+4) + j*nzz + k*nxx*nzz]) +
                       1029.0f*(Txz[(i+3) + j*nzz + k*nxx*nzz] - Txz[(i-2) + j*nzz + k*nxx*nzz]) +
                       8575.0f*(Txz[(i-1) + j*nzz + k*nxx*nzz] - Txz[(i+2) + j*nzz + k*nxx*nzz]) +
                     128625.0f*(Txz[(i+1) + j*nzz + k*nxx*nzz] - Txz[i + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

        float Bx = 0.5f*(B[i + (j+1)*nzz + k*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);

        Vx[index] += dt*Bx*(dTxx_dx + dTxy_dy + dTxz_dz); 
    }

    if((i >= 3) && (i < nzz-3) && (j >= 3) && (j < nxx-4) && (k > 3) && (k < nyy-3)) 
    {
        float dTxy_dx = (75.0f*(Txy[i + (j-3)*nzz + k*nxx*nzz] - Txy[i + (j+4)*nzz + k*nxx*nzz]) +
                       1029.0f*(Txy[i + (j+3)*nzz + k*nxx*nzz] - Txy[i + (j-2)*nzz + k*nxx*nzz]) +
                       8575.0f*(Txy[i + (j-1)*nzz + k*nxx*nzz] - Txy[i + (j+2)*nzz + k*nxx*nzz]) +
                     128625.0f*(Txy[i + (j+1)*nzz + k*nxx*nzz] - Txy[i + j*nzz + k*nxx*nzz])) / (107520.0f*dx);

        float dTyy_dy = (75.0f*(Tyy[i + j*nzz + (k-4)*nxx*nzz] - Tyy[i + j*nzz + (k+3)*nxx*nzz]) +
                       1029.0f*(Tyy[i + j*nzz + (k+2)*nxx*nzz] - Tyy[i + j*nzz + (k-3)*nxx*nzz]) +
                       8575.0f*(Tyy[i + j*nzz + (k-2)*nxx*nzz] - Tyy[i + j*nzz + (k+1)*nxx*nzz]) +
                     128625.0f*(Tyy[i + j*nzz + k*nxx*nzz]     - Tyy[i + j*nzz + (k-1)*nxx*nzz])) / (107520.0f*dy);

        float dTyz_dz = (75.0f*(Tyz[(i-3) + j*nzz + k*nxx*nzz] - Tyz[(i+4) + j*nzz + k*nxx*nzz]) +
                       1029.0f*(Tyz[(i+3) + j*nzz + k*nxx*nzz] - Tyz[(i-2) + j*nzz + k*nxx*nzz]) +
                       8575.0f*(Tyz[(i-1) + j*nzz + k*nxx*nzz] - Tyz[(i+2) + j*nzz + k*nxx*nzz]) +
                     128625.0f*(Tyz[(i+1) + j*nzz + k*nxx*nzz] - Tyz[i + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

        float By = 0.5f*(B[i + j*nzz + (k+1)*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);

        Vy[index] += dt*By*(dTxy_dx + dTyy_dy + dTyz_dz); 
    }    

    if((i > 3) && (i < nzz-3) && (j >= 3) && (j < nxx-4) && (k >= 3) && (k < nyy-4)) 
    {
        float dTxz_dx = (75.0f*(Txz[i + (j-3)*nzz + k*nxx*nzz] - Txz[i + (j+4)*nzz + k*nxx*nzz]) +
                       1029.0f*(Txz[i + (j+3)*nzz + k*nxx*nzz] - Txz[i + (j-2)*nzz + k*nxx*nzz]) +
                       8575.0f*(Txz[i + (j-1)*nzz + k*nxx*nzz] - Txz[i + (j+2)*nzz + k*nxx*nzz]) +
                     128625.0f*(Txz[i + (j+1)*nzz + k*nxx*nzz] - Txz[i + j*nzz + k*nxx*nzz])) / (107520.0f*dx);

        float dTyz_dy = (75.0f*(Tyz[i + j*nzz + (k-3)*nxx*nzz] - Tyz[i + j*nzz + (k+4)*nxx*nzz]) +
                       1029.0f*(Tyz[i + j*nzz + (k+3)*nxx*nzz] - Tyz[i + j*nzz + (k-2)*nxx*nzz]) +
                       8575.0f*(Tyz[i + j*nzz + (k-1)*nxx*nzz] - Tyz[i + j*nzz + (k+2)*nxx*nzz]) +
                     128625.0f*(Tyz[i + j*nzz + (k+1)*nxx*nzz] - Tyz[i + j*nzz + k*nxx*nzz])) / (107520.0f*dy);

        float dTzz_dz = (75.0f*(Tzz[(i-4) + j*nzz + k*nxx*nzz] - Tzz[(i+3) + j*nzz + k*nxx*nzz]) +
                       1029.0f*(Tzz[(i+2) + j*nzz + k*nxx*nzz] - Tzz[(i-3) + j*nzz + k*nxx*nzz]) +
                       8575.0f*(Tzz[(i-2) + j*nzz + k*nxx*nzz] - Tzz[(i+1) + j*nzz + k*nxx*nzz]) +
                     128625.0f*(Tzz[i + j*nzz + k*nxx*nzz]     - Tzz[(i-1) + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

        float Bz = 0.5f*(B[(i+1) + j*nzz + k*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);

        Vz[index] += dt*Bz*(dTxz_dx + dTyz_dy + dTzz_dz); 
    }

	float damper = get_boundary_damper(damp1D, damp2D, damp3D, i, j, k, nxx, nyy, nzz, nb);

    if (index < nxx*nyy*nzz) 
    {
        Vx[index] *= damper;
        Vy[index] *= damper;
        Vz[index] *= damper;

        Txx[index] *= damper;
        Tyy[index] *= damper;
        Tzz[index] *= damper;
        Txz[index] *= damper;
        Tyz[index] *= damper;
        Txy[index] *= damper;    
    }
}

__global__ void compute_seismogram(float * P, int * rIdx, int * rIdy, int * rIdz, float * seismogram, int spread, int tId, int tlag, int nt, int nxx, int nzz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if ((index < spread) && (tId >= tlag))
        seismogram[(tId - tlag) + index*nt] = P[rIdz[index] + rIdx[index]*nzz + rIdy[index]*nxx*nzz];
}

__device__ float get_boundary_damper(float * damp1D, float * damp2D, float * damp3D, int i, int j, int k, int nxx, int nyy, int nzz, int nabc)
{
    float damper;

    // global case
    if ((i >= nabc) && (i < nzz-nabc) && (j >= nabc) && (j < nxx-nabc) && (k >= nabc) && (k < nyy-nabc))
    {
        damper = 1.0f;
    }

    // 1D damping
    else if((i < nabc) && (j >= nabc) && (j < nxx-nabc) && (k >= nabc) && (k < nyy-nabc)) 
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

    return damper;
}
