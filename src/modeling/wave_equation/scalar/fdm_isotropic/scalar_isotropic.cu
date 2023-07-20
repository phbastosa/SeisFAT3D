# include "scalar_isotropic.cuh"

void Scalar_Isotropic::set_model_parameters()
{
    modeling_method = std::string("scalar_isotropic");
    modeling_message = std::string("[3] - SConstant density acoustic isotropic media\n\n");

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
        
        compute_pressure<<<blocksPerGrid,threadsPerBlock>>>(Pressure,U_pre,U_pas,dtVp2,damp1D,damp2D,damp3D,wavelet,sId,time_id,dx,dy,dz,nxx,nyy,nzz,nb,nbzu);
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

    cudaFree(damp1D);
    cudaFree(damp2D);
    cudaFree(damp3D);
    
    cudaFree(wavelet);
    cudaFree(Pressure);
    cudaFree(seismogram);
}

__global__ void compute_pressure(float * Pressure, float * U_pre, float * U_pas, float * dtVp2, float * damp1D, float * damp2D, float * damp3D, float * wavelet, int sId, int time_id, float dx, float dy, float dz, int nxx, int nyy, int nzz, int nb, int nbzu)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         // y direction
    int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
    int i = (int) (index - j*nzz - k*nxx*nzz); // z direction
    
    if (index == 0) U_pre[sId] += wavelet[time_id] / (dx*dy*dz);

    if((i >= 4) && (i < nzz-4) && (j >= 4) && (j < nxx-4) && (k >= 4) && (k < nyy-4)) 
    {
        float d2_Px2 = (- 9.0f*(U_pre[i + (j-4)*nzz + k*nxx*nzz] + U_pre[i + (j+4)*nzz + k*nxx*nzz])
                    +   128.0f*(U_pre[i + (j-3)*nzz + k*nxx*nzz] + U_pre[i + (j+3)*nzz + k*nxx*nzz])
                    -  1008.0f*(U_pre[i + (j-2)*nzz + k*nxx*nzz] + U_pre[i + (j+2)*nzz + k*nxx*nzz])
                    +  8064.0f*(U_pre[i + (j-1)*nzz + k*nxx*nzz] + U_pre[i + (j+1)*nzz + k*nxx*nzz])
                    - 14350.0f*(U_pre[i + j*nzz + k*nxx*nzz]))/(5040.0f*powf(dx, 2.0f));

        float d2_Py2 = (- 9.0f*(U_pre[i + j*nzz + (k-4)*nxx*nzz] + U_pre[i + j*nzz + (k+4)*nxx*nzz])
                    +   128.0f*(U_pre[i + j*nzz + (k-3)*nxx*nzz] + U_pre[i + j*nzz + (k+3)*nxx*nzz])
                    -  1008.0f*(U_pre[i + j*nzz + (k-2)*nxx*nzz] + U_pre[i + j*nzz + (k+2)*nxx*nzz])
                    +  8064.0f*(U_pre[i + j*nzz + (k-1)*nxx*nzz] + U_pre[i + j*nzz + (k+1)*nxx*nzz])
                    - 14350.0f*(U_pre[i + j*nzz + k*nxx*nzz]))/(5040.0f*powf(dy,2.0f));

        float d2_Pz2 = (- 9.0f*(U_pre[(i-4) + j*nzz + k*nxx*nzz] + U_pre[(i+4) + j*nzz + k*nxx*nzz])
                    +   128.0f*(U_pre[(i-3) + j*nzz + k*nxx*nzz] + U_pre[(i+3) + j*nzz + k*nxx*nzz])
                    -  1008.0f*(U_pre[(i-2) + j*nzz + k*nxx*nzz] + U_pre[(i+2) + j*nzz + k*nxx*nzz])
                    +  8064.0f*(U_pre[(i-1) + j*nzz + k*nxx*nzz] + U_pre[(i+1) + j*nzz + k*nxx*nzz])
                    - 14350.0f*(U_pre[i + j*nzz + k*nxx*nzz]))/(5040.0f*powf(dz,2.0f));
    
        
        Pressure[index] = dtVp2[index] * (d2_Px2 + d2_Py2 + d2_Pz2) + 2.0f*U_pre[index] - U_pas[index]; 
    }

    float damper = get_boundary_damper(damp1D,damp2D,damp3D,i,j,k,nxx,nyy,nzz,nb,nbzu);
    
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