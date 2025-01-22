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

    T = new float[volsize]();
    TT = new float[nPoints]();

    current_xrec = new int[max_spread]();
    current_yrec = new int[max_spread]();
    current_zrec = new int[max_spread]();

    synthetic_data = new float[nt*max_spread]();

    cudaMalloc((void**)&(d_P), volsize*sizeof(float));
    cudaMalloc((void**)&(d_T), volsize*sizeof(float));

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

    float t0 = 2.0f*sqrtf(M_PI) / fmax;
    float fc = fmax / (3.0f * sqrtf(M_PI));

    tlag = (int)(t0 / dt) + 1;

    for (int n = 0; n < nt; n++)
    {
        float td = n*dt - t0;

        float arg = M_PI*M_PI*M_PI*fc*fc*td*td;

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

void Elastic::forward_solver()
{
    eikonal->srcId = srcId;
    eikonal->forward_solver();

    eikonal->reduce_boundary(eikonal->T, TT);
    
    expand_boundary(TT, T);
    
    cudaMemcpy(d_T, T, volsize*sizeof(float), cudaMemcpyHostToDevice);

    initialization();

    propagation();

    cudaMemcpy(synthetic_data, seismogram, nt*geometry->spread[srcId]*sizeof(float), cudaMemcpyDeviceToHost);
}

void Elastic::export_synthetic_data()
{
    std::string data_file = data_folder + modeling_type + "_nStations" + std::to_string(geometry->spread[srcId]) + "_nSamples" + std::to_string(nt) + "_shot_" + std::to_string(geometry->sInd[srcId]+1) + ".bin";
    export_binary_float(data_file, synthetic_data, nt*geometry->spread[srcId]);    
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

__global__ void compute_velocity(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * B, float * T,  float * damp1D, float * damp2D, float * damp3D, float dx, float dy, float dz, float dt, int tId, int tlag, int nxx, int nyy, int nzz, int nb)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         
    int j = (int) (index - k*nxx*nzz) / nzz;   
    int i = (int) (index - j*nzz - k*nxx*nzz); 

    if ((index < nxx*nyy*nzz) && (T[index] < (float)(tId + tlag)*dt))
    {
        if((i >= 3) && (i < nzz-4) && (j > 3) && (j < nxx-3) && (k >= 3) && (k < nyy-4)) 
        {
            float dTxx_dx = (FDM1*(Txx[i + (j-4)*nzz + k*nxx*nzz] - Txx[i + (j+3)*nzz + k*nxx*nzz]) +
                             FDM2*(Txx[i + (j+2)*nzz + k*nxx*nzz] - Txx[i + (j-3)*nzz + k*nxx*nzz]) +
                             FDM3*(Txx[i + (j-2)*nzz + k*nxx*nzz] - Txx[i + (j+1)*nzz + k*nxx*nzz]) +
                             FDM4*(Txx[i + j*nzz + k*nxx*nzz]     - Txx[i + (j-1)*nzz + k*nxx*nzz])) / dx;

            float dTxy_dy = (FDM1*(Txy[i + j*nzz + (k-3)*nxx*nzz] - Txy[i + j*nzz + (k+4)*nxx*nzz]) +
                             FDM2*(Txy[i + j*nzz + (k+3)*nxx*nzz] - Txy[i + j*nzz + (k-2)*nxx*nzz]) +
                             FDM3*(Txy[i + j*nzz + (k-1)*nxx*nzz] - Txy[i + j*nzz + (k+2)*nxx*nzz]) +
                             FDM4*(Txy[i + j*nzz + (k+1)*nxx*nzz] - Txy[i + j*nzz + k*nxx*nzz])) / dy;

            float dTxz_dz = (FDM1*(Txz[(i-3) + j*nzz + k*nxx*nzz] - Txz[(i+4) + j*nzz + k*nxx*nzz]) +
                             FDM2*(Txz[(i+3) + j*nzz + k*nxx*nzz] - Txz[(i-2) + j*nzz + k*nxx*nzz]) +
                             FDM3*(Txz[(i-1) + j*nzz + k*nxx*nzz] - Txz[(i+2) + j*nzz + k*nxx*nzz]) +
                             FDM4*(Txz[(i+1) + j*nzz + k*nxx*nzz] - Txz[i + j*nzz + k*nxx*nzz])) / dz;

            float Bx = 0.5f*(B[i + (j+1)*nzz + k*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);

            Vx[index] += dt*Bx*(dTxx_dx + dTxy_dy + dTxz_dz); 
        }

        if((i >= 3) && (i < nzz-3) && (j >= 3) && (j < nxx-4) && (k > 3) && (k < nyy-3)) 
        {
            float dTxy_dx = (FDM1*(Txy[i + (j-3)*nzz + k*nxx*nzz] - Txy[i + (j+4)*nzz + k*nxx*nzz]) +
                             FDM2*(Txy[i + (j+3)*nzz + k*nxx*nzz] - Txy[i + (j-2)*nzz + k*nxx*nzz]) +
                             FDM3*(Txy[i + (j-1)*nzz + k*nxx*nzz] - Txy[i + (j+2)*nzz + k*nxx*nzz]) +
                             FDM4*(Txy[i + (j+1)*nzz + k*nxx*nzz] - Txy[i + j*nzz + k*nxx*nzz])) / dx;

            float dTyy_dy = (FDM1*(Tyy[i + j*nzz + (k-4)*nxx*nzz] - Tyy[i + j*nzz + (k+3)*nxx*nzz]) +
                             FDM2*(Tyy[i + j*nzz + (k+2)*nxx*nzz] - Tyy[i + j*nzz + (k-3)*nxx*nzz]) +
                             FDM3*(Tyy[i + j*nzz + (k-2)*nxx*nzz] - Tyy[i + j*nzz + (k+1)*nxx*nzz]) +
                             FDM4*(Tyy[i + j*nzz + k*nxx*nzz]     - Tyy[i + j*nzz + (k-1)*nxx*nzz])) / dy;

            float dTyz_dz = (FDM1*(Tyz[(i-3) + j*nzz + k*nxx*nzz] - Tyz[(i+4) + j*nzz + k*nxx*nzz]) +
                             FDM2*(Tyz[(i+3) + j*nzz + k*nxx*nzz] - Tyz[(i-2) + j*nzz + k*nxx*nzz]) +
                             FDM3*(Tyz[(i-1) + j*nzz + k*nxx*nzz] - Tyz[(i+2) + j*nzz + k*nxx*nzz]) +
                             FDM4*(Tyz[(i+1) + j*nzz + k*nxx*nzz] - Tyz[i + j*nzz + k*nxx*nzz])) / dz;

            float By = 0.5f*(B[i + j*nzz + (k+1)*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);

            Vy[index] += dt*By*(dTxy_dx + dTyy_dy + dTyz_dz); 
        }    

        if((i > 3) && (i < nzz-3) && (j >= 3) && (j < nxx-4) && (k >= 3) && (k < nyy-4)) 
        {
            float dTxz_dx = (FDM1*(Txz[i + (j-3)*nzz + k*nxx*nzz] - Txz[i + (j+4)*nzz + k*nxx*nzz]) +
                             FDM2*(Txz[i + (j+3)*nzz + k*nxx*nzz] - Txz[i + (j-2)*nzz + k*nxx*nzz]) +
                             FDM3*(Txz[i + (j-1)*nzz + k*nxx*nzz] - Txz[i + (j+2)*nzz + k*nxx*nzz]) +
                             FDM4*(Txz[i + (j+1)*nzz + k*nxx*nzz] - Txz[i + j*nzz + k*nxx*nzz])) / dx;

            float dTyz_dy = (FDM1*(Tyz[i + j*nzz + (k-3)*nxx*nzz] - Tyz[i + j*nzz + (k+4)*nxx*nzz]) +
                             FDM2*(Tyz[i + j*nzz + (k+3)*nxx*nzz] - Tyz[i + j*nzz + (k-2)*nxx*nzz]) +
                             FDM3*(Tyz[i + j*nzz + (k-1)*nxx*nzz] - Tyz[i + j*nzz + (k+2)*nxx*nzz]) +
                             FDM4*(Tyz[i + j*nzz + (k+1)*nxx*nzz] - Tyz[i + j*nzz + k*nxx*nzz])) / dy;

            float dTzz_dz = (FDM1*(Tzz[(i-4) + j*nzz + k*nxx*nzz] - Tzz[(i+3) + j*nzz + k*nxx*nzz]) +
                             FDM2*(Tzz[(i+2) + j*nzz + k*nxx*nzz] - Tzz[(i-3) + j*nzz + k*nxx*nzz]) +
                             FDM3*(Tzz[(i-2) + j*nzz + k*nxx*nzz] - Tzz[(i+1) + j*nzz + k*nxx*nzz]) +
                             FDM4*(Tzz[i + j*nzz + k*nxx*nzz]     - Tzz[(i-1) + j*nzz + k*nxx*nzz])) / dz;

            float Bz = 0.5f*(B[(i+1) + j*nzz + k*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);

            Vz[index] += dt*Bz*(dTxz_dx + dTyz_dy + dTzz_dz); 
        }

    	float damper = get_boundary_damper(damp1D, damp2D, damp3D, i, j, k, nxx, nyy, nzz, nb);

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
