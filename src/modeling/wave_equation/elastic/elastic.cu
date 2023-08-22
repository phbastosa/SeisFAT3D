# include "elastic.cuh"

void Elastic::set_parameters()
{
    general_modeling_parameters();
    wave_modeling_parameters();

    set_acquisition_geometry();

    set_velocity_model();
    set_density_model();
    set_shear_model();

    set_boundaries();
    set_model_boundaries();

    set_gridded_geometry();

    set_modeling_volumes();

    set_wavelet();
    set_dampers();
    set_outputs();
}

void Elastic::set_density_model()
{
    Rho = new float[nPoints]();

    import_binary_float(catch_parameter("rho_model_file", file), Rho, nPoints);
}

void Elastic::set_shear_model()
{
    Vs = new float[nPoints]();

    import_binary_float(catch_parameter("vs_model_file", file), Vs, nPoints);
}

void Elastic::set_model_boundaries()
{
    float * b = new float[volsize]();
    float * m = new float[volsize]();
    float * l = new float[volsize]();

    expand_boundary(Rho, b);
    expand_boundary(Vs, m);
    expand_boundary(Vp, l);

    for (int index = 0; index < volsize; index++)
    {
        m[index] = b[index]*powf(m[index], 2.0f);
        l[index] = b[index]*powf(l[index], 2.0f) - 2.0f*m[index];
        b[index] = 1.0f / b[index];
    }
    
	cudaMalloc((void**)&(B), volsize*sizeof(float));
	cudaMalloc((void**)&(M), volsize*sizeof(float));
	cudaMalloc((void**)&(L), volsize*sizeof(float));

	cudaMemcpy(B, b, volsize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(M, m, volsize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(L, l, volsize*sizeof(float), cudaMemcpyHostToDevice);    

    delete[] b;
    delete[] m;
    delete[] l;
}

void Elastic::set_modeling_volumes()
{
    modeling_method = std::string("elastic");
    modeling_message = std::string("[5] - Elastic isotropic media\n\n");

	cudaMalloc((void**)&(Vx), volsize*sizeof(float));
	cudaMalloc((void**)&(Vy), volsize*sizeof(float));
	cudaMalloc((void**)&(Vz), volsize*sizeof(float));

	cudaMalloc((void**)&(Txx), volsize*sizeof(float));
	cudaMalloc((void**)&(Tyy), volsize*sizeof(float));
	cudaMalloc((void**)&(Tzz), volsize*sizeof(float));
	cudaMalloc((void**)&(Txz), volsize*sizeof(float));
	cudaMalloc((void**)&(Tyz), volsize*sizeof(float));
	cudaMalloc((void**)&(Txy), volsize*sizeof(float));

    cudaMalloc((void**)&(Pressure), volsize*sizeof(float));
}

void Elastic::set_wavelet()
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

    if (import_wavelet) 
        import_binary_float(wavelet_file, aux, nt);

	cudaMalloc((void**)&(wavelet), nt*sizeof(float));
	cudaMemcpy(wavelet, aux, nt*sizeof(float), cudaMemcpyHostToDevice);

    delete[] aux;
    delete[] ricker;
}

void Elastic::info_message()
{
    general_modeling_message();
    
    set_modeling_message();
}

void Elastic::initial_setup()
{
    cudaMemset(Vx, 0.0f, volsize*sizeof(float));
    cudaMemset(Vy, 0.0f, volsize*sizeof(float));
    cudaMemset(Vz, 0.0f, volsize*sizeof(float));

    cudaMemset(Txx, 0.0f, volsize*sizeof(float));
    cudaMemset(Tyy, 0.0f, volsize*sizeof(float));
    cudaMemset(Tzz, 0.0f, volsize*sizeof(float));
    cudaMemset(Txz, 0.0f, volsize*sizeof(float));
    cudaMemset(Tyz, 0.0f, volsize*sizeof(float));
    cudaMemset(Txy, 0.0f, volsize*sizeof(float));

    cudaMemset(Pressure, 0.0f, volsize*sizeof(float));
    cudaMemset(seismogram, 0.0f, nt*total_nodes*sizeof(float));

    int sidx = (int)(geometry->shots.x[shot_id] / dx) + nbxl;
    int sidy = (int)(geometry->shots.y[shot_id] / dy) + nbyl;
    int sidz = (int)(geometry->shots.z[shot_id] / dz) + nbzu;

    source_id = sidz + sidx*nzz + sidy*nxx*nzz;

    isnap = 0;
}

void Elastic::forward_solver()
{
    for (time_id = 0; time_id < nt; time_id++)
    {
        show_progress();
        
        compute_velocity<<<blocksPerGrid,threadsPerBlock>>>(Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,B,wavelet,source_id,time_id,dx,dy,dz,dt,nxx,nyy,nzz);
        cudaDeviceSynchronize();

        compute_stress<<<blocksPerGrid,threadsPerBlock>>>(Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,Pressure,M,L,damp1D,damp2D,damp3D,dx,dy,dz,dt,nxx,nyy,nzz,nb,nbzu);
        cudaDeviceSynchronize();
    
        get_snapshots();
        get_seismogram();
    }
}

void Elastic::free_space()
{
    cudaFree(B);
    cudaFree(M);
    cudaFree(L);

    cudaFree(Vx);
    cudaFree(Vy);
    cudaFree(Vz);

    cudaFree(Txx);
    cudaFree(Tyy);
    cudaFree(Tzz);
    cudaFree(Txz);
    cudaFree(Tyz);
    cudaFree(Txy);

    cudaFree(damp1D);
    cudaFree(damp2D);
    cudaFree(damp3D);
    
    cudaFree(wavelet);
    cudaFree(Pressure);
    cudaFree(seismogram);
}

__global__ void compute_velocity(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * B, float * wavelet, int source_id, int time_id, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         // y direction
    int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
    int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

    if (index == 0)
    {
        Txx[source_id] += wavelet[time_id] / (dx*dy*dz);
        Tyy[source_id] += wavelet[time_id] / (dx*dy*dz);
        Tzz[source_id] += wavelet[time_id] / (dx*dy*dz);
    }

    if((i >= 3) && (i < nzz-4) && (j > 3) && (j < nxx-3) && (k >= 3) && (k < nyy-4)) 
    {
        float d_Txx_dx = (75.0f*(Txx[i + (j-4)*nzz + k*nxx*nzz] - Txx[i + (j+3)*nzz + k*nxx*nzz]) +
                        1029.0f*(Txx[i + (j+2)*nzz + k*nxx*nzz] - Txx[i + (j-3)*nzz + k*nxx*nzz]) +
                        8575.0f*(Txx[i + (j-2)*nzz + k*nxx*nzz] - Txx[i + (j+1)*nzz + k*nxx*nzz]) +
                      128625.0f*(Txx[i + j*nzz + k*nxx*nzz]     - Txx[i + (j-1)*nzz + k*nxx*nzz])) / (107520.0f*dx);

        float d_Txy_dy = (75.0f*(Txy[i + j*nzz + (k-3)*nxx*nzz] - Txy[i + j*nzz + (k+4)*nxx*nzz]) +
                        1029.0f*(Txy[i + j*nzz + (k+3)*nxx*nzz] - Txy[i + j*nzz + (k-2)*nxx*nzz]) +
                        8575.0f*(Txy[i + j*nzz + (k-1)*nxx*nzz] - Txy[i + j*nzz + (k+2)*nxx*nzz]) +
                      128625.0f*(Txy[i + j*nzz + (k+1)*nxx*nzz] - Txy[i + j*nzz + k*nxx*nzz])) / (107520.0f*dy);

        float d_Txz_dz = (75.0f*(Txz[(i-3) + j*nzz + k*nxx*nzz] - Txz[(i+4) + j*nzz + k*nxx*nzz]) +
                        1029.0f*(Txz[(i+3) + j*nzz + k*nxx*nzz] - Txz[(i-2) + j*nzz + k*nxx*nzz]) +
                        8575.0f*(Txz[(i-1) + j*nzz + k*nxx*nzz] - Txz[(i+2) + j*nzz + k*nxx*nzz]) +
                      128625.0f*(Txz[(i+1) + j*nzz + k*nxx*nzz] - Txz[i + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

        float Bx = 0.5f*(B[i + (j+1)*nzz + k*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);

        Vx[index] += dt*Bx*(d_Txx_dx + d_Txy_dy + d_Txz_dz); 
    }

    if((i >= 3) && (i < nzz-3) && (j >= 3) && (j < nxx-4) && (k > 3) && (k < nyy-3)) 
    {
        float d_Txy_dx = (75.0f*(Txy[i + (j-3)*nzz + k*nxx*nzz] - Txy[i + (j+4)*nzz + k*nxx*nzz]) +
                        1029.0f*(Txy[i + (j+3)*nzz + k*nxx*nzz] - Txy[i + (j-2)*nzz + k*nxx*nzz]) +
                        8575.0f*(Txy[i + (j-1)*nzz + k*nxx*nzz] - Txy[i + (j+2)*nzz + k*nxx*nzz]) +
                      128625.0f*(Txy[i + (j+1)*nzz + k*nxx*nzz] - Txy[i + j*nzz + k*nxx*nzz])) / (107520.0f*dx);

        float  d_Tyy_dy = (75.0f*(Tyy[i + j*nzz + (k-4)*nxx*nzz] - Tyy[i + j*nzz + (k+3)*nxx*nzz]) +
                         1029.0f*(Tyy[i + j*nzz + (k+2)*nxx*nzz] - Tyy[i + j*nzz + (k-3)*nxx*nzz]) +
                         8575.0f*(Tyy[i + j*nzz + (k-2)*nxx*nzz] - Tyy[i + j*nzz + (k+1)*nxx*nzz]) +
                       128625.0f*(Tyy[i + j*nzz + k*nxx*nzz]     - Tyy[i + j*nzz + (k-1)*nxx*nzz])) / (107520.0f*dy);

        float d_Tyz_dz = (75.0f*(Tyz[(i-3) + j*nzz + k*nxx*nzz] - Tyz[(i+4) + j*nzz + k*nxx*nzz]) +
                        1029.0f*(Tyz[(i+3) + j*nzz + k*nxx*nzz] - Tyz[(i-2) + j*nzz + k*nxx*nzz]) +
                        8575.0f*(Tyz[(i-1) + j*nzz + k*nxx*nzz] - Tyz[(i+2) + j*nzz + k*nxx*nzz]) +
                      128625.0f*(Tyz[(i+1) + j*nzz + k*nxx*nzz] - Tyz[i + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

        float By = 0.5f*(B[i + j*nzz + (k+1)*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);

        Vy[index] += dt*By*(d_Txy_dx + d_Tyy_dy + d_Tyz_dz); 
    }    

    if((i > 3) && (i < nzz-3) && (j >= 3) && (j < nxx-4) && (k >= 3) && (k < nyy-4)) 
    {
        float d_Txz_dx = (75.0f*(Txz[i + (j-3)*nzz + k*nxx*nzz] - Txz[i + (j+4)*nzz + k*nxx*nzz]) +
                        1029.0f*(Txz[i + (j+3)*nzz + k*nxx*nzz] - Txz[i + (j-2)*nzz + k*nxx*nzz]) +
                        8575.0f*(Txz[i + (j-1)*nzz + k*nxx*nzz] - Txz[i + (j+2)*nzz + k*nxx*nzz]) +
                      128625.0f*(Txz[i + (j+1)*nzz + k*nxx*nzz] - Txz[i + j*nzz + k*nxx*nzz])) / (107520.0f*dx);

        float d_Tyz_dy = (75.0f*(Tyz[i + j*nzz + (k-3)*nxx*nzz] - Tyz[i + j*nzz + (k+4)*nxx*nzz]) +
                        1029.0f*(Tyz[i + j*nzz + (k+3)*nxx*nzz] - Tyz[i + j*nzz + (k-2)*nxx*nzz]) +
                        8575.0f*(Tyz[i + j*nzz + (k-1)*nxx*nzz] - Tyz[i + j*nzz + (k+2)*nxx*nzz]) +
                      128625.0f*(Tyz[i + j*nzz + (k+1)*nxx*nzz] - Tyz[i + j*nzz + k*nxx*nzz])) / (107520.0f*dy);

        float d_Tzz_dz = (75.0f*(Tzz[(i-4) + j*nzz + k*nxx*nzz] - Tzz[(i+3) + j*nzz + k*nxx*nzz]) +
                        1029.0f*(Tzz[(i+2) + j*nzz + k*nxx*nzz] - Tzz[(i-3) + j*nzz + k*nxx*nzz]) +
                        8575.0f*(Tzz[(i-2) + j*nzz + k*nxx*nzz] - Tzz[(i+1) + j*nzz + k*nxx*nzz]) +
                      128625.0f*(Tzz[i + j*nzz + k*nxx*nzz]     - Tzz[(i-1) + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

        float Bz = 0.5f*(B[(i+1) + j*nzz + k*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);

        Vz[index] += dt*Bz*(d_Txz_dx + d_Tyz_dy + d_Tzz_dz); 
    }
}

__global__ void compute_stress(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * Pressure, float * M, float * L, float * damp1D, float * damp2D, float * damp3D, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz, int nb, int nbzu)
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

        Txx[index] += dt*((L[index] + 2*M[index])*d_Vx_dx + L[index]*(d_Vy_dy + d_Vz_dz));
        Tyy[index] += dt*((L[index] + 2*M[index])*d_Vy_dy + L[index]*(d_Vx_dx + d_Vz_dz));
        Tzz[index] += dt*((L[index] + 2*M[index])*d_Vz_dz + L[index]*(d_Vx_dx + d_Vy_dy));                    
    }

    if((i >= 0) && (i < nzz) && (j > 3) && (j < nxx-3) && (k > 3) && (k < nyy-3)) 
    {
        float d_Vx_dy = (75.0f*(Vx[i + j*nzz + (k-4)*nxx*nzz] - Vx[i + j*nzz + (k+3)*nxx*nzz]) +
                       1029.0f*(Vx[i + j*nzz + (k+2)*nxx*nzz] - Vx[i + j*nzz + (k-3)*nxx*nzz]) +
                       8575.0f*(Vx[i + j*nzz + (k-2)*nxx*nzz] - Vx[i + j*nzz + (k+1)*nxx*nzz]) +
                     128625.0f*(Vx[i + j*nzz + k*nxx*nzz]     - Vx[i + j*nzz + (k-1)*nxx*nzz])) / (107520.0f*dy);

        float d_Vy_dx = (75.0f*(Vy[i + (j-4)*nzz + k*nxx*nzz] - Vy[i + (j+3)*nzz + k*nxx*nzz]) +
                       1029.0f*(Vy[i + (j+2)*nzz + k*nxx*nzz] - Vy[i + (j-3)*nzz + k*nxx*nzz]) +
                       8575.0f*(Vy[i + (j-2)*nzz + k*nxx*nzz] - Vy[i + (j+1)*nzz + k*nxx*nzz]) +
                     128625.0f*(Vy[i + j*nzz + k*nxx*nzz]     - Vy[i + (j-1)*nzz + k*nxx*nzz])) / (107520.0f*dx);

        float M_xy = powf(0.25*(1/M[i + (j+1)*nzz + (k+1)*nxx*nzz] + 1/M[i + (j+1)*nzz + k*nxx*nzz] + 
                                1/M[i + j*nzz + (k+1)*nxx*nzz]     + 1/M[i + j*nzz + k*nxx*nzz]), -1.0f);

        Txy[index] += dt*M_xy*(d_Vx_dy + d_Vy_dx);
    }

    if((i > 3) && (i < nzz-3) && (j > 3) && (j < nxx-3) && (k >= 0) && (k < nyy)) 
    {
        float d_Vx_dz = (75.0f*(Vx[(i-4) + j*nzz + k*nxx*nzz] - Vx[(i+3) + j*nzz + k*nxx*nzz]) +
                       1029.0f*(Vx[(i+2) + j*nzz + k*nxx*nzz] - Vx[(i-3) + j*nzz + k*nxx*nzz]) +
                       8575.0f*(Vx[(i-2) + j*nzz + k*nxx*nzz] - Vx[(i+1) + j*nzz + k*nxx*nzz]) +
                     128625.0f*(Vx[i + j*nzz + k*nxx*nzz]     - Vx[(i-1) + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

        float d_Vz_dx = (75.0f*(Vz[i + (j-4)*nzz + k*nxx*nzz] - Vz[i + (j+3)*nzz + k*nxx*nzz]) +
                       1029.0f*(Vz[i + (j+2)*nzz + k*nxx*nzz] - Vz[i + (j-3)*nzz + k*nxx*nzz]) +
                       8575.0f*(Vz[i + (j-2)*nzz + k*nxx*nzz] - Vz[i + (j+1)*nzz + k*nxx*nzz]) +
                     128625.0f*(Vz[i + j*nzz + k*nxx*nzz]     - Vz[i + (j-1)*nzz + k*nxx*nzz])) / (107520.0f*dx);

        float M_xz = powf(0.25*(1/M[(i+1) + (j+1)*nzz + k*nxx*nzz] + 1/M[i + (j+1)*nzz + k*nxx*nzz] + 
                                1/M[(i+1) + j*nzz + k*nxx*nzz]     + 1/M[i + j*nzz + k*nxx*nzz]), -1.0f);

        Txz[index] += dt*M_xz*(d_Vx_dz + d_Vz_dx);
    }

    if((i > 3) && (i < nzz-3) && (j >= 0) && (j < nxx) && (k > 3) && (k < nyy-3)) 
    {
        float d_Vy_dz = (75.0f*(Vy[(i-4) + j*nzz + k*nxx*nzz] - Vy[(i+3) + j*nzz + k*nxx*nzz]) +
                       1029.0f*(Vy[(i+2) + j*nzz + k*nxx*nzz] - Vy[(i-3) + j*nzz + k*nxx*nzz]) +
                       8575.0f*(Vy[(i-2) + j*nzz + k*nxx*nzz] - Vy[(i+1) + j*nzz + k*nxx*nzz]) +
                     128625.0f*(Vy[i + j*nzz + k*nxx*nzz]     - Vy[(i-1) + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

        float d_Vz_dy = (75.0f*(Vz[i + j*nzz + (k-4)*nxx*nzz] - Vz[i + j*nzz + (k+3)*nxx*nzz]) +
                       1029.0f*(Vz[i + j*nzz + (k+2)*nxx*nzz] - Vz[i + j*nzz + (k-3)*nxx*nzz]) +
                       8575.0f*(Vz[i + j*nzz + (k-2)*nxx*nzz] - Vz[i + j*nzz + (k+1)*nxx*nzz]) +
                     128625.0f*(Vz[i + j*nzz + k*nxx*nzz]     - Vz[i + j*nzz + (k-1)*nxx*nzz])) / (107520.0f*dy);

        float M_yz = powf(0.25*(1/M[(i+1) + j*nzz + (k+1)*nxx*nzz] + 1/M[i + j*nzz + (k+1)*nxx*nzz] + 
                                1/M[(i+1) + j*nzz + k*nxx*nzz] +     1/M[i + j*nzz + k*nxx*nzz]), -1.0f);

        Tyz[index] += dt*M_yz*(d_Vy_dz + d_Vz_dy);
    }
        
    float damper = get_boundary_damper(damp1D,damp2D,damp3D,i,j,k,nxx,nyy,nzz,nb,nbzu);

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

        Pressure[index] = (Txx[index] + Tyy[index] + Tzz[index]) / 3.0f;
    }
}
