# include "acoustic.cuh"

void Acoustic::set_parameters()
{
    general_modeling_parameters();
    wave_modeling_parameters();

    set_acquisition_geometry();
    set_gridded_geometry();

    set_velocity_model();
    set_density_model();

    set_boundaries();
    set_model_boundaries();

    set_modeling_volumes();

    set_wavelet();
    set_dampers();
    set_outputs();
}

void Acoustic::set_density_model()
{
    Rho = new float[nPoints]();

    // import_binary_float(catch_parameter("rho_model_file", file), Rho, nPoints);

    for (int index = 0; index < nPoints; index++) Rho[index] = 1000.0f; 
}

void Acoustic::set_model_boundaries()
{
    float * b = new float[volsize]();
    float * k = new float[volsize]();

    expand_boundary(Vp, k);
    expand_boundary(Rho, b);

    for (int index = 0; index < volsize; index++)
    {
        k[index] = k[index]*k[index]*b[index]; 
        b[index] = 1.0f / b[index]; 
    }

	cudaMalloc((void**)&(B), volsize*sizeof(float));
	cudaMalloc((void**)&(K), volsize*sizeof(float));

	cudaMemcpy(B, b, volsize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(K, k, volsize*sizeof(float), cudaMemcpyHostToDevice);

    delete[] b;
    delete[] k;
}

void Acoustic::set_modeling_volumes()
{
    modeling_method = std::string("acoustic");
    modeling_message = std::string("[4] - Variable density acoustic isotropic media\n\n");

	cudaMalloc((void**)&(Vx), volsize*sizeof(float));
	cudaMalloc((void**)&(Vy), volsize*sizeof(float));
	cudaMalloc((void**)&(Vz), volsize*sizeof(float));

    cudaMalloc((void**)&(Pressure), volsize*sizeof(float));
}

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

    if (import_wavelet) import_binary_float(wavelet_file, aux, nt);

	cudaMalloc((void**)&(wavelet), nt*sizeof(float));
	cudaMemcpy(wavelet, aux, nt*sizeof(float), cudaMemcpyHostToDevice);

    delete[] aux;
    delete[] ricker;
}

void Acoustic::info_message()
{
    general_modeling_message();
    
    set_modeling_message();
}

void Acoustic::initial_setup()
{
    cudaMemset(Vx, 0.0f, volsize*sizeof(float));
    cudaMemset(Vy, 0.0f, volsize*sizeof(float));
    cudaMemset(Vz, 0.0f, volsize*sizeof(float));

    cudaMemset(Pressure, 0.0f, volsize*sizeof(float));
    cudaMemset(seismogram, 0.0f, nt*total_nodes*sizeof(float));

    int sidx = (int)(geometry->shots.x[shot_id] / dx) + nbxl;
    int sidy = (int)(geometry->shots.y[shot_id] / dy) + nbyl;
    int sidz = (int)(geometry->shots.z[shot_id] / dz) + nbzu;

    source_id = sidz + sidx*nzz + sidy*nxx*nzz;

    isnap = 0;
}

void Acoustic::forward_solver()
{
    for (time_id = 0; time_id < nt; time_id++)
    {
        show_progress();
        
        compute_velocity<<<blocksPerGrid,threadsPerBlock>>>(Pressure,Vx,Vy,Vz,B,wavelet,source_id,time_id,dx,dy,dz,dt,nxx,nyy,nzz);
        cudaDeviceSynchronize();

        compute_pressure<<<blocksPerGrid,threadsPerBlock>>>(Pressure,Vx,Vy,Vz,K,damp1D,damp2D,damp3D,dx,dy,dz,dt,nxx,nyy,nzz,nbzu,nb);
        cudaDeviceSynchronize();
        
        get_snapshots();
        get_seismogram();
    }
}

void Acoustic::free_space()
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

__global__ void compute_velocity(float * Pressure, float * Vx, float * Vy, float * Vz, float * B, float * wavelet, int source_id, int time_id, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         // y direction
    int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
    int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

    if (index == 0) Pressure[source_id] += wavelet[time_id] / (dx*dy*dz);

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












