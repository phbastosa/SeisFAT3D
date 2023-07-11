# include "acoustic_isotropic.cuh"

void Acoustic_Isotropic::set_model_parameters()
{
    modeling_method = std::string("acoustic_isotropic");
    modeling_message = std::string("[4] - Solution for the wave equation in variable density acoustic isotropic media\n\n");

    float * vp = new float[nPoints]();
    float * rho = new float[nPoints]();

    // import_binary_float(catch_parameter("vp_model_file", file), vp, nPoints);
    // import_binary_float(catch_parameter("rho_model_file", file), rho, nPoints);

    for (int index = 0; index < nPoints; index++) 
    {
        vp[index] = 1500.0f;
        rho[index] = 1000.0f;
    }   

    Vp = new float[volsize]();
    Rho = new float[volsize]();

    expand_boundary(vp, Vp);
    expand_boundary(rho, Rho);

    delete[] vp;
    delete[] rho;

    float * b = new float[volsize]();
    float * k = new float[volsize]();

    for (int index = 0; index < volsize; index++)
    {
        b[index] = 1.0f / Rho[index];
        k[index] = Vp[index]*Vp[index]*Rho[index];
    }
    
	cudaMalloc((void**)&(B), volsize*sizeof(float));
	cudaMalloc((void**)&(K), volsize*sizeof(float));

	cudaMemcpy(B, b, volsize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(K, k, volsize*sizeof(float), cudaMemcpyHostToDevice);

    delete[] b;
    delete[] k;
}

void Acoustic_Isotropic::set_wavefields()
{
	cudaMalloc((void**)&(Vx), volsize*sizeof(float));
	cudaMalloc((void**)&(Vy), volsize*sizeof(float));
	cudaMalloc((void**)&(Vz), volsize*sizeof(float));

    cudaMalloc((void**)&(Pressure), volsize*sizeof(float));
}

void Acoustic_Isotropic::initial_setup()
{
    cudaMemset(Vx, 0.0f, volsize*sizeof(float));
    cudaMemset(Vy, 0.0f, volsize*sizeof(float));
    cudaMemset(Vz, 0.0f, volsize*sizeof(float));

    cudaMemset(Pressure, 0.0f, volsize*sizeof(float));
    cudaMemset(seismogram, 0.0f, nt*total_nodes*sizeof(float));

    int sidx = (int)(geometry->shots.x[shot_id] / dx) + nbxl;
    int sidy = (int)(geometry->shots.y[shot_id] / dy) + nbyl;
    int sidz = (int)(geometry->shots.z[shot_id] / dz) + nbzu;

    sId = sidz + sidx*nzz + sidy*nxx*nzz;

    isnap = 0;
}

void Acoustic_Isotropic::forward_solver()
{
    for (time_id = 0; time_id < nt; time_id++)
    {
        show_progress();
        
        compute_velocity<<<blocksPerGrid,threadsPerBlock>>>(Pressure,Vx,Vy,Vz,B,wavelet,sId,time_id,dx,dy,dz,dt,nxx,nyy,nzz);
        cudaDeviceSynchronize();

        compute_pressure<<<blocksPerGrid,threadsPerBlock>>>(Pressure,Vx,Vy,Vz,K,damp1D,damp2D,damp3D,dx,dy,dz,dt,nxx,nyy,nzz,nbzu,nb);
        cudaDeviceSynchronize();
        
        get_snapshots();
        get_seismogram();
    }
}

void Acoustic_Isotropic::free_space()
{
    cudaFree(B);
    cudaFree(K);

    cudaFree(Vx);
    cudaFree(Vy);
    cudaFree(Vz);

    cudaFree(damp1D);
    cudaFree(damp2D);
    cudaFree(damp3D);
    
    cudaFree(wavelet);
    cudaFree(Pressure);    
    cudaFree(seismogram);
}

__global__ void compute_velocity(float * Pressure, float * Vx, float * Vy, float * Vz, float * B, float * wavelet, int sId, int time_id, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         // y direction
    int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
    int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

    if (index == 0) Pressure[sId] += wavelet[time_id] / (dx*dy*dz);

    if((i > 3) && (i < nzz-3) && (j > 3) && (j < nxx-3) && (k > 3) && (k < nyy-3)) 
    {
        float d_P_dx = (75.0f*(Pressure[i + (j-4)*nzz + k*nxx*nzz] - Pressure[i + (j+3)*nzz + k*nxx*nzz]) +
                      1029.0f*(Pressure[i + (j+2)*nzz + k*nxx*nzz] - Pressure[i + (j-3)*nzz + k*nxx*nzz]) +
                      8575.0f*(Pressure[i + (j-2)*nzz + k*nxx*nzz] - Pressure[i + (j+1)*nzz + k*nxx*nzz]) +
                    128625.0f*(Pressure[i + j*nzz + k*nxx*nzz]     - Pressure[i + (j-1)*nzz + k*nxx*nzz])) / (107520.0f*dx);

        float d_P_dy = (75.0f*(Pressure[i + j*nzz + (k-4)*nxx*nzz] - Pressure[i + j*nzz + (k+3)*nxx*nzz]) +
                      1029.0f*(Pressure[i + j*nzz + (k+2)*nxx*nzz] - Pressure[i + j*nzz + (k-3)*nxx*nzz]) +
                      8575.0f*(Pressure[i + j*nzz + (k-2)*nxx*nzz] - Pressure[i + j*nzz + (k+1)*nxx*nzz]) +
                    128625.0f*(Pressure[i + j*nzz + k*nxx*nzz]     - Pressure[i + j*nzz + (k-1)*nxx*nzz])) / (107520.0f*dy);

        float d_P_dz = (75.0f*(Pressure[(i-4) + j*nzz + k*nxx*nzz] - Pressure[(i+3) + j*nzz + k*nxx*nzz]) +
                      1029.0f*(Pressure[(i+2) + j*nzz + k*nxx*nzz] - Pressure[(i-3) + j*nzz + k*nxx*nzz]) +
                      8575.0f*(Pressure[(i-2) + j*nzz + k*nxx*nzz] - Pressure[(i+1) + j*nzz + k*nxx*nzz]) +
                    128625.0f*(Pressure[i + j*nzz + k*nxx*nzz]     - Pressure[(i-1) + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

        float Bx = 0.5f*(B[i + (j+1)*nzz + k*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);
        float By = 0.5f*(B[i + j*nzz + (k+1)*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);
        float Bz = 0.5f*(B[(i+1) + j*nzz + k*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);

        Vx[index] += -dt*Bx*d_P_dx; 
        Vy[index] += -dt*By*d_P_dy; 
        Vz[index] += -dt*Bz*d_P_dz; 
    }
}

__global__ void compute_pressure(float * Pressure, float * Vx, float * Vy, float * Vz, float * K, float * damp1D, float * damp2D, float * damp3D, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz, int nbzu, int nb)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         // y direction
    int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
    int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

    if((i >= 3) && (i < nzz-4) && (j >= 3) && (j < nxx-4) && (k >= 3) && (k < nyy-4)) 
    {
        float d_Vx_dx = (75.0f*(Vx[i + (j-3)*nzz + k*nxx*nzz] - Vx[i + (j+4)*nzz + k*nxx*nzz]) +
                       1029.0f*(Vx[i + (j+3)*nzz + k*nxx*nzz] - Vx[i + (j-2)*nzz + k*nxx*nzz]) +
                       8575.0f*(Vx[i + (j-1)*nzz + k*nxx*nzz] - Vx[i + (j+2)*nzz + k*nxx*nzz]) +
                     128625.0f*(Vx[i + (j+1)*nzz + k*nxx*nzz] - Vx[i + j*nzz + k*nxx*nzz])) / (107520.0f*dx);

        float d_Vy_dy = (75.0f*(Vy[i + j*nzz + (k-3)*nxx*nzz] - Vy[i + j*nzz + (k+4)*nxx*nzz]) +
                       1029.0f*(Vy[i + j*nzz + (k+3)*nxx*nzz] - Vy[i + j*nzz + (k-2)*nxx*nzz]) +
                       8575.0f*(Vy[i + j*nzz + (k-1)*nxx*nzz] - Vy[i + j*nzz + (k+2)*nxx*nzz]) +
                     128625.0f*(Vy[i + j*nzz + (k+1)*nxx*nzz] - Vy[i + j*nzz + k*nxx*nzz])) / (107520.0f*dy);

        float d_Vz_dz = (75.0f*(Vz[(i-3) + j*nzz + k*nxx*nzz] - Vz[(i+4) + j*nzz + k*nxx*nzz]) +
                       1029.0f*(Vz[(i+3) + j*nzz + k*nxx*nzz] - Vz[(i-2) + j*nzz + k*nxx*nzz]) +
                       8575.0f*(Vz[(i-1) + j*nzz + k*nxx*nzz] - Vz[(i+2) + j*nzz + k*nxx*nzz]) +
                     128625.0f*(Vz[(i+1) + j*nzz + k*nxx*nzz] - Vz[i + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

        Pressure[index] += -dt*K[index]*(d_Vx_dx + d_Vy_dy + d_Vz_dz); 
    }

    float damper = get_boundary_damper(damp1D,damp2D,damp3D,i,j,k,nxx,nyy,nzz,nb,nbzu);

    if (index < nxx*nyy*nzz)
    {
        Vx[index] *= damper;
        Vy[index] *= damper;
        Vz[index] *= damper;
        
        Pressure[index] *= damper;
    }   
}