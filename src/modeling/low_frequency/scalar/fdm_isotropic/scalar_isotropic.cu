# include "scalar_isotropic.cuh"

void Scalar_Isotropic::set_modeling_message()
{
    std::cout<<"Running:\n";
    std::cout<<"[3] - Solution for wave equation in constant density acoustic isotropic media\n"; 
    std::cout<<"    Cerjan et al. (1985)\n\n";

    std::cout<<"Modeling progress: " << floorf(100.0f * (float)(time_id+1) / (float)(nt)) <<" %\n\n";
}

void Scalar_Isotropic::set_model_parameters()
{
    modeling_method = std::string("scalar_isotropic");

    float * vp = new float[nPoints]();
    float * dtvp2 = new float[volsize]();

    // import_binary_float(catch_parameter("vp_model_file", file), vp, nPoints);

    for (int index = 0; index < nPoints; index++) 
    {
        vp[index] = 1500.0f;
    }

    Vp = new float[volsize]();

    expand_boundary(vp, Vp);

    delete[] vp;

    for (int index = 0; index < volsize; index++) 
    {
        dtvp2[index] = dt*dt*Vp[index]*Vp[index];
    }

	cudaMalloc((void**)&(dtVp2), volsize*sizeof(float));

	cudaMemcpy(dtVp2, dtvp2, volsize*sizeof(float), cudaMemcpyHostToDevice);
    
    delete[] dtvp2;
}

void Scalar_Isotropic::set_wavefields()
{
	cudaMalloc((void**)&(U_pre), volsize*sizeof(float));
	cudaMalloc((void**)&(U_pas), volsize*sizeof(float));
    cudaMalloc((void**)&(Pressure), volsize*sizeof(float));
}

void Scalar_Isotropic::initial_setup()
{
    cudaMemset(U_pre, 0.0f, volsize*sizeof(float));
    cudaMemset(U_pas, 0.0f, volsize*sizeof(float));
    cudaMemset(Pressure, 0.0f, volsize*sizeof(float));
    cudaMemset(seismogram, 0.0f, nt*total_nodes*sizeof(float));

    int sidx = (int)(geometry->shots.x[shot_id] / dx) + nbxl;
    int sidy = (int)(geometry->shots.y[shot_id] / dy) + nbyl;
    int sidz = (int)(geometry->shots.z[shot_id] / dz) + nbzu;

    sId = sidz + sidx*nzz + sidy*nxx*nzz;

    isnap = 0;
}

void Scalar_Isotropic::forward_solver()
{
    for (time_id = 0; time_id < nt; time_id++)
    {
        show_progress();
        
        apply_wavelet<<<1,1>>>(U_pre,wavelet,sId,time_id,dx,dy,dz);
        cudaDeviceSynchronize();

        compute_pressure<<<blocksPerGrid,threadsPerBlock>>>(Pressure,U_pre,U_pas,dtVp2,damp1D,damp2D,damp3D,dx,dy,dz,nxx,nyy,nzz,nb,nbzu);
        cudaDeviceSynchronize();

        update_pressure<<<blocksPerGrid,threadsPerBlock>>>(Pressure,U_pre,U_pas,volsize);
        cudaDeviceSynchronize();
    
        get_snapshots();
        get_seismogram();
    }
}

void Scalar_Isotropic::free_space()
{
    cudaFree(dtVp2);
    cudaFree(U_pre);
    cudaFree(U_pas);

    cudaFree(Pressure);

    cudaFree(damp1D);
    cudaFree(damp2D);
    cudaFree(damp3D);
    
    cudaFree(wavelet);
    
    cudaFree(seismogram);
}

__global__ void apply_wavelet(float * Pressure, float * wavelet, int sId, int time_id, float dx, float dy, float dz)
{
    Pressure[sId] += wavelet[time_id] / (dx*dy*dz);
}

__global__ void compute_pressure(float * Pressure, float * U_pre, float * U_pas, float * dtVp2, float * damp1D, float * damp2D, float * damp3D, float dx, float dy, float dz, int nxx, int nyy, int nzz, int nb, int nbzu)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         // y direction
    int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
    int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

    float damper = 1.0f;

    float d2_Px2, d2_Py2, d2_Pz2;

    if((i >= 4) && (i < nzz-4) && (j >= 4) && (j < nxx-4) && (k >= 4) && (k < nyy-4)) 
    {
        d2_Px2 = (- 9.0f*(U_pre[i + (j-4)*nzz + k*nxx*nzz] + U_pre[i + (j+4)*nzz + k*nxx*nzz])
              +   128.0f*(U_pre[i + (j-3)*nzz + k*nxx*nzz] + U_pre[i + (j+3)*nzz + k*nxx*nzz])
              -  1008.0f*(U_pre[i + (j-2)*nzz + k*nxx*nzz] + U_pre[i + (j+2)*nzz + k*nxx*nzz])
              +  8064.0f*(U_pre[i + (j-1)*nzz + k*nxx*nzz] + U_pre[i + (j+1)*nzz + k*nxx*nzz])
              - 14350.0f*(U_pre[i + j*nzz + k*nxx*nzz]))/(5040.0f*powf(dx, 2.0f));

        d2_Py2 = (- 9.0f*(U_pre[i + j*nzz + (k-4)*nxx*nzz] + U_pre[i + j*nzz + (k+4)*nxx*nzz])
              +   128.0f*(U_pre[i + j*nzz + (k-3)*nxx*nzz] + U_pre[i + j*nzz + (k+3)*nxx*nzz])
              -  1008.0f*(U_pre[i + j*nzz + (k-2)*nxx*nzz] + U_pre[i + j*nzz + (k+2)*nxx*nzz])
              +  8064.0f*(U_pre[i + j*nzz + (k-1)*nxx*nzz] + U_pre[i + j*nzz + (k+1)*nxx*nzz])
              - 14350.0f*(U_pre[i + j*nzz + k*nxx*nzz]))/(5040.0f*powf(dy,2.0f));

        d2_Pz2 = (- 9.0f*(U_pre[(i-4) + j*nzz + k*nxx*nzz] + U_pre[(i+4) + j*nzz + k*nxx*nzz])
              +   128.0f*(U_pre[(i-3) + j*nzz + k*nxx*nzz] + U_pre[(i+3) + j*nzz + k*nxx*nzz])
              -  1008.0f*(U_pre[(i-2) + j*nzz + k*nxx*nzz] + U_pre[(i+2) + j*nzz + k*nxx*nzz])
              +  8064.0f*(U_pre[(i-1) + j*nzz + k*nxx*nzz] + U_pre[(i+1) + j*nzz + k*nxx*nzz])
              - 14350.0f*(U_pre[i + j*nzz + k*nxx*nzz]))/(5040.0f*powf(dz,2.0f));
    
        
        Pressure[index] = dtVp2[index] * (d2_Px2 + d2_Py2 + d2_Pz2) + 2.0f*U_pre[index] - U_pas[index]; 
    }
        
    // 1D damping
    if((i < nbzu) && (j >= nb) && (j < nxx-nb) && (k >= nb) && (k < nyy-nb)) 
    {
        damper = damp1D[i];
    }         
    else if((i >= nzz-nb) && (i < nzz) && (j >= nb) && (j < nxx-nb) && (k >= nb) && (k < nyy-nb)) 
    {
        damper = damp1D[nb-(i-(nzz-nb))-1];
    }         
    else if((i >= nbzu) && (i < nzz-nb) && (j >= 0) && (j < nb) && (k >= nb) && (k < nyy-nb)) 
    {
        damper = damp1D[j];
    }
    else if((i >= nbzu) && (i < nzz-nb) && (j >= nxx-nb) && (j < nxx) && (k >= nb) && (k < nyy-nb)) 
    {
        damper = damp1D[nb-(j-(nxx-nb))-1];
    }
    else if((i >= nbzu) && (i < nzz-nb) && (j >= nb) && (j < nxx-nb) && (k >= 0) && (k < nb)) 
    {
        damper = damp1D[k];
    }
    else if((i >= nbzu) && (i < nzz-nb) && (j >= nb) && (j < nxx-nb) && (k >= nyy-nb) && (k < nyy)) 
    {
        damper = damp1D[nb-(k-(nyy-nb))-1];
    }

    // 2D damping 
    else if((i >= nbzu) && (i < nzz-nb) && (j >= 0) && (j < nb) && (k >= 0) && (k < nb))
    {
        damper = damp2D[j + k*nb];
    }
    else if((i >= nbzu) && (i < nzz-nb) && (j >= nxx-nb) && (j < nxx) && (k >= 0) && (k < nb))
    {
        damper = damp2D[nb-(j-(nxx-nb))-1 + k*nb];
    }
    else if((i >= nbzu) && (i < nzz-nb) && (j >= 0) && (j < nb) && (k >= nyy-nb) && (k < nyy))
    {
        damper = damp2D[j + (nb-(k-(nyy-nb))-1)*nb];
    }
    else if((i >= nbzu) && (i < nzz-nb) && (j >= nxx-nb) && (j < nxx) && (k >= nyy-nb) && (k < nyy))
    {
        damper = damp2D[nb-(j-(nxx-nb))-1 + (nb-(k-(nyy-nb))-1)*nb];
    }

    else if((i >= 0) && (i < nb) && (j >= nb) && (j < nxx-nb) && (k >= 0) && (k < nb))
    {
        damper = damp2D[i + k*nb];
    }
    else if((i >= nzz-nb) && (i < nzz) && (j >= nb) && (j < nxx-nb) && (k >= 0) && (k < nb))
    {
        damper = damp2D[nb-(i-(nzz-nb))-1 + k*nb];
    }
    else if((i >= 0) && (i < nb) && (j >= nb) && (j < nxx-nb) && (k >= nyy-nb) && (k < nyy))
    {
        damper = damp2D[i + (nb-(k-(nyy-nb))-1)*nb];
    }
    else if((i >= nzz-nb) && (i < nzz) && (j >= nb) && (j < nxx-nb) && (k >= nyy-nb) && (k < nyy))
    {
        damper = damp2D[nb-(i-(nzz-nb))-1 + (nb-(k-(nyy-nb))-1)*nb];
    }

    else if((i >= 0) && (i < nb) && (j >= 0) && (j < nb) && (k >= nb) && (k < nyy-nb))
    {
        damper = damp2D[i + j*nb];
    }
    else if((i >= nzz-nb) && (i < nzz) && (j >= 0) && (j < nb) && (k >= nb) && (k < nyy-nb))
    {
        damper = damp2D[nb-(i-(nzz-nb))-1 + j*nb];
    }
    else if((i >= 0) && (i < nb) && (j >= nxx-nb) && (j < nxx) && (k >= nb) && (k < nyy-nb))
    {
        damper = damp2D[i + (nb-(j-(nxx-nb))-1)*nb];
    }
    else if((i >= nzz-nb) && (i < nzz) && (j >= nxx-nb) && (j < nxx) && (k >= nb) && (k < nyy-nb))
    {
        damper = damp2D[nb-(i-(nzz-nb))-1 + (nb-(j-(nxx-nb))-1)*nb];
    }

    // 3D damping
    else if((i >= 0) && (i < nb) && (j >= 0) && (j < nb) && (k >= 0) && (k < nb))
    {
        damper = damp3D[i + j*nb + k*nb*nb];
    }
    else if((i >= nzz-nb) && (i < nzz) && (j >= 0) && (j < nb) && (k >= 0) && (k < nb))
    {
        damper = damp3D[nb-(i-(nzz-nb))-1 + j*nb + k*nb*nb];
    }
    else if((i >= 0) && (i < nb) && (j >= nxx-nb) && (j < nxx) && (k >= 0) && (k < nb))
    {
        damper = damp3D[i + (nb-(j-(nxx-nb))-1)*nb + k*nb*nb];
    }
    else if((i >= 0) && (i < nb) && (j >= 0) && (j < nb) && (k >= nyy-nb) && (k < nyy))
    {
        damper = damp3D[i + j*nb + (nb-(k-(nyy-nb))-1)*nb*nb];
    }
    else if((i >= nzz-nb) && (i < nzz) && (j >= nxx-nb) && (j < nxx) && (k >= 0) && (k < nb))
    {
        damper = damp3D[nb-(i-(nzz-nb))-1 + (nb-(j-(nxx-nb))-1)*nb + k*nb*nb];
    }
    else if((i >= nzz-nb) && (i < nzz) && (j >= 0) && (j < nb) && (k >= nyy-nb) && (k < nyy))
    {
        damper = damp3D[nb-(i-(nzz-nb))-1 + j*nb + (nb-(k-(nyy-nb))-1)*nb*nb];
    }
    else if((i >= 0) && (i < nb) && (j >= nxx-nb) && (j < nxx) && (k >= nyy-nb) && (k < nyy))
    {
        damper = damp3D[i + (nb-(j-(nxx-nb))-1)*nb + (nb-(k-(nyy-nb))-1)*nb*nb];
    }
    else if((i >= nzz-nb) && (i < nzz) && (j >= nxx-nb) && (j < nxx) && (k >= nyy-nb) && (k < nyy))
    {
        damper = damp3D[nb-(i-(nzz-nb))-1 + (nb-(j-(nxx-nb))-1)*nb + (nb-(k-(nyy-nb))-1)*nb*nb];
    }
    
    if (index < nxx*nyy*nzz)
    {
        U_pas[index] *= damper;    
        U_pre[index] *= damper;
        Pressure[index] *= damper;
    }
}

__global__ void update_pressure(float * Pressure, float * U_pre, float * U_pas, int volsize)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (index < volsize)
    {
        U_pas[index] = U_pre[index];        
        U_pre[index] = Pressure[index];
    }
}