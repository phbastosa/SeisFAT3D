# include "elastic.cuh"

void Elastic::set_model_parameters()
{
    nxx = nx + 2*nb;
    nyy = ny + 2*nb;
    nzz = nz + 2*nb;

    nPoints = nx*ny*nz;
    volsize = nxx*nyy*nzz;

    totalBlocks = (int)(volsize / threadsPerBlock);

    Vp = new float[nPoints]();
    Vs = new float[nPoints]();
    Rho = new float[nPoints]();

    h_B = new float[volsize]();
    h_M = new float[volsize]();
    h_L = new float[volsize]();

    import_binary_float(catch_parameter("vp_model_file", file), Vp, nPoints);
    import_binary_float(catch_parameter("vs_model_file", file), Vs, nPoints);
    import_binary_float(catch_parameter("rho_model_file", file), Rho, nPoints);

    int total_time = std::stoi(catch_parameter("total_time", file));

    dt = std::stof(catch_parameter("time_spacing", file));

    nt = (int)(ceilf(total_time / dt)) + 1; 

    set_wavelet();
}

void Elastic::set_wavelet()
{
    float * wavelet = new float[nt]();

    float fmax = std::stof(catch_parameter("max_frequency", file));
    float tlag = std::stof(catch_parameter("wavelet_shift", file));
    float amp = std::stof(catch_parameter("max_amplitude", file));

    float pi = 4.0f * atanf(1.0f);  

    float fc = fmax / (3.0f * sqrtf(pi));

    float * aux = new float[nt]();

    for (int n = 0; n < nt; n++)
    {        
        float arg = pi*((n*dt - tlag)*fc*pi)*((n*dt - tlag)*fc*pi);
        
        aux[n] = (1 - 2*arg) * expf(-arg);    

        float sum = 0.0f;
        for (int k = 0; k < n+1; k++)     
            sum += aux[k];
    
        wavelet[n] = amp * sum;
    }

	cudaMalloc((void**)&(d_wavelet), nt*sizeof(float));
	cudaMemcpy(d_wavelet, wavelet, nt*sizeof(float), cudaMemcpyHostToDevice);

    delete[] aux;
    delete[] wavelet;
}

void Elastic::set_model_boundaries()
{
    for (int z = nb; z < nzz - nb; z++)
    {
        for (int y = nb; y < nyy - nb; y++)
        {
            for (int x = nb; x < nxx - nb; x++)
            {
                int ind1 = (z - nb) + (x - nb)*nz + (y - nb)*nx*nz;
                int ind2 = z + x*nzz + y*nxx*nzz;

                h_B[ind2] = 1.0f / Rho[ind1];
                h_M[ind2] = Rho[ind1]*powf(Vs[ind1], 2.0f);
                h_L[ind2] = Rho[ind1]*powf(Vp[ind1], 2.0f) - 2.0f*h_M[ind2]; 
            }
        }
    }

    for (int z = 0; z < nb; z++)
    {
        for (int y = nb; y < nyy - nb; y++)
        {
            for (int x = nb; x < nxx - nb; x++)
            {
                int ind1 = z + x*nzz + y*nxx*nzz;
                int ind2 = (nzz - z - 1) + x*nzz + y*nxx*nzz;
                int ind3 = 0 + (x - nb)*nz + (y - nb)*nx*nz;
                int ind4 = (nz - 1) + (x - nb)*nz + (y - nb)*nx*nz;         
  
                h_B[ind1] = 1.0f / Rho[ind3];
                h_B[ind2] = 1.0f / Rho[ind4];

                h_M[ind1] = Rho[ind3]*powf(Vs[ind3], 2.0f);
                h_M[ind2] = Rho[ind4]*powf(Vs[ind4], 2.0f);

                h_L[ind1] = Rho[ind3]*powf(Vp[ind3], 2.0f) - 2.0f*h_M[ind1]; 
                h_L[ind2] = Rho[ind4]*powf(Vp[ind4], 2.0f) - 2.0f*h_M[ind2]; 
            }
        }
    }

    for (int x = 0; x < nb; x++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int y = nb; y < nyy - nb; y++)
            {
                h_B[z + x*nzz + y*nxx*nzz] = h_B[z + nb*nzz + y*nxx*nzz];
                h_B[z + (nxx - x - 1)*nzz + y*nxx*nzz] = h_B[z + (nxx - nb - 1)*nzz + y*nxx*nzz];

                h_M[z + x*nzz + y*nxx*nzz] = h_M[z + nb*nzz + y*nxx*nzz];
                h_M[z + (nxx - x - 1)*nzz + y*nxx*nzz] = h_M[z + (nxx - nb - 1)*nzz + y*nxx*nzz];

                h_L[z + x*nzz + y*nxx*nzz] = h_L[z + nb*nzz + y*nxx*nzz];
                h_L[z + (nxx - x - 1)*nzz + y*nxx*nzz] = h_L[z + (nxx - nb - 1)*nzz + y*nxx*nzz];
            }
        }
    }

    for (int y = 0; y < nb; y++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int x = 0; x < nxx; x++)
            {
                h_B[z + x*nzz + y*nxx*nzz] = h_B[z + x*nzz + nb*nxx*nzz];
                h_B[z + x*nzz + (nyy - y - 1)*nxx*nzz] = h_B[z + x*nzz + (nyy - nb - 1)*nxx*nzz];

                h_M[z + x*nzz + y*nxx*nzz] = h_M[z + x*nzz + nb*nxx*nzz];
                h_M[z + x*nzz + (nyy - y - 1)*nxx*nzz] = h_M[z + x*nzz + (nyy - nb - 1)*nxx*nzz];

                h_L[z + x*nzz + y*nxx*nzz] = h_L[z + x*nzz + nb*nxx*nzz];
                h_L[z + x*nzz + (nyy - y - 1)*nxx*nzz] = h_L[z + x*nzz + (nyy - nb - 1)*nxx*nzz];
            }
        }
    }

	cudaMalloc((void**)&(B), volsize*sizeof(float));
	cudaMalloc((void**)&(M), volsize*sizeof(float));
	cudaMalloc((void**)&(L), volsize*sizeof(float));

	cudaMemcpy(B, h_B, volsize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(M, h_M, volsize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(L, h_L, volsize*sizeof(float), cudaMemcpyHostToDevice);

    delete[] h_B;
    delete[] h_M;
    delete[] h_L;


    set_dampers();
}

void Elastic::set_dampers()
{
    float * damp1D = new float[nb]();
    float * damp2D = new float[nb*nb]();
    float * damp3D = new float[nb*nb*nb]();

    float factor = std::stoi(catch_parameter("boundary_damper", file));

    for (int i = 0; i < nb; i++) 
    {
        damp1D[i] = expf(-powf(factor * (nb - i), 2.0f));
    }

    for(int i = 0; i < nb; i++) 
    {
        for (int j = 0; j < nb; j++)
        {   
            damp2D[j + i*nb] += damp1D[i]; // up to bottom
            damp2D[i + j*nb] += damp1D[i]; // left to right
        }
    }

    for (int i  = 0; i < nb; i++)
    {
        for(int j = 0; j < nb; j++)
        {
            for(int k = 0; k < nb; k++)
            {
                damp3D[i + j*nb + k*nb*nb] += damp2D[i + j*nb]; // XY plane
                damp3D[i + j*nb + k*nb*nb] += damp2D[j + k*nb]; // ZX plane
                damp3D[i + j*nb + k*nb*nb] += damp2D[i + k*nb]; // ZY plane
            }
        }
    }    

    for (int index = 0; index < nb*nb; index++)
        damp2D[index] -= 1.0f;

    for (int index = 0; index < nb*nb*nb; index++)
        damp3D[index] -= 5.0f;    

	cudaMalloc((void**)&(d_damp1D), nb*sizeof(float));
	cudaMalloc((void**)&(d_damp2D), nb*nb*sizeof(float));
	cudaMalloc((void**)&(d_damp3D), nb*nb*nb*sizeof(float));

	cudaMemcpy(d_damp1D, damp1D, nb*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_damp2D, damp2D, nb*nb*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_damp3D, damp3D, nb*nb*nb*sizeof(float), cudaMemcpyHostToDevice);

    delete[] damp1D;
    delete[] damp2D;
    delete[] damp3D;
}

void Elastic::set_wavefields()
{
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

void Elastic::set_outputs()
{
    nsnap = std::stoi(catch_parameter("n_snapshots", file));
    dsnap = (int)(nt / (nsnap-1));

    receiver_output_samples = nt*total_nodes;
    wavefield_output_samples = nPoints*nsnap;

    if (export_receiver_output)
    {
        cudaMalloc((void**)&(seismogram), receiver_output_samples*sizeof(float));
        receiver_output = new float[receiver_output_samples]();
    }

    if (export_wavefield_output)
    {   
        snapshot = new float[volsize]();
        wavefield_output = new float[wavefield_output_samples]();
    }
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

    int sidx = (int)(geometry->shots.x[shot_id] / dx) + nb;
    int sidy = (int)(geometry->shots.y[shot_id] / dy) + nb;
    int sidz = (int)(geometry->shots.z[shot_id] / dz) + nb;

    sId = sidz + sidx*nzz + sidy*nxx*nzz;

    isnap = 0;
}

void Elastic::forward_solver()
{
    for (time_id = 0; time_id < nt; time_id++)
    {
        show_progress();
        
        apply_wavelet<<<1,1>>>(Txx, Tyy, Tzz, d_wavelet, sId, time_id, dx, dy, dz);
        cudaDeviceSynchronize();

        compute_velocity<<<totalBlocks,threadsPerBlock>>>(Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,B,dx,dy,dz,dt,nxx,nyy,nzz);
        cudaDeviceSynchronize();

        compute_stress<<<totalBlocks,threadsPerBlock>>>(Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,M,L,d_damp1D,d_damp2D,d_damp3D,dx,dy,dz,dt,nxx,nyy,nzz,nb);
        cudaDeviceSynchronize();

        compute_pressure<<<totalBlocks,threadsPerBlock>>>(Txx,Tyy,Tzz,Pressure,volsize);
        cudaDeviceSynchronize();

        // // if (export_receiver_output)
        // //     get_seismogram
    
        if (export_wavefield_output)
        {
            if (nt % dsnap == 0)
            {
                cudaMemcpy(snapshot, Pressure, volsize*sizeof(float), cudaMemcpyDeviceToHost);       

                # pragma omp parallel for 
                for (int index = 0; index < nPoints; index++)
                {
                    int y = (int) (index / (nx*nz));         
                    int x = (int) (index - y*nx*nz) / nz;    
                    int z = (int) (index - x*nz - y*nx*nz);  

                    wavefield_output[z + x*nz + y*nx*nz + isnap*nPoints] = snapshot[(z + nb) + (x + nb)*nzz + (y + nb)*nxx*nzz];
                }

                isnap++;
            }
        }
    }
}

void Elastic::show_progress()
{
    if (time_id % (nt / 100) == 0)        
    {
        Modeling::info_message();

        std::cout<<"Modeling progress: " << floorf(100.0f * (float)(time_id+1) / (float)(nt)) <<" %"<<std::endl;
    }
}

void Elastic::build_outputs()
{
    cudaMemcpy(seismogram, receiver_output, receiver_output_samples*sizeof(float), cudaMemcpyDeviceToHost);

    receiver_output_file = receiver_output_folder + "elastic_pressure_" + std::to_string(nt) + "x" + std::to_string(geometry->nodes.total) + "_shot_" + std::to_string(shot_id+1) + ".bin";

    wavefield_output_file = wavefield_output_folder + "elastic_pressure_" + std::to_string(nz) + "x" + std::to_string(nx) + "x" + std::to_string(ny) + "_" + std::to_string(nsnap) + "_"+ "_shot_" + std::to_string(shot_id+1) + ".bin";
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

    cudaFree(Pressure);
    cudaFree(d_damp1D);
    cudaFree(d_damp2D);
    cudaFree(d_damp3D);
    cudaFree(d_wavelet);
    cudaFree(seismogram);
}

__global__ void apply_wavelet(float * Txx, float * Tyy, float * Tzz, float * wavelet, int sId, int time_id, float dx, float dy, float dz)
{    
    Txx[sId] += wavelet[time_id] / (dx*dy*dz);
    Tyy[sId] += wavelet[time_id] / (dx*dy*dz);
    Tzz[sId] += wavelet[time_id] / (dx*dy*dz);
}

__global__ void compute_velocity(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * B, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         // y direction
    int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
    int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

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
                      128625.0f*(Tzz[i + j*nzz + k*nxx*nzz]     - Tzz[(i+1) + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

        float Bz = 0.5f*(B[(i+1) + j*nzz + k*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);

        Vz[index] += dt*Bz*(d_Txz_dx + d_Tyz_dy + d_Tzz_dz); 
    }
}

__global__ void compute_stress(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * M, float * L, float * damp1D, float * damp2D, float * damp3D, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz, int nb)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         // y direction
    int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
    int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

    float damper = 1.0f;

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
        
    // 1D damping
    if((i < nb) && (j >= nb) && (j < nxx-nb) && (k >= nb) && (k < nyy-nb)) 
    {
        damper = damp1D[i];
    }         
    else if((i >= nzz-nb) && (i < nzz) && (j >= nb) && (j < nxx-nb) && (k >= nb) && (k < nyy-nb)) 
    {
        damper = damp1D[nb-(i-(nzz-nb))-1];
    }         
    else if((i >= nb) && (i < nzz-nb) && (j >= 0) && (j < nb) && (k >= nb) && (k < nyy-nb)) 
    {
        damper = damp1D[j];
    }
    else if((i >= nb) && (i < nzz-nb) && (j >= nxx-nb) && (j < nxx) && (k >= nb) && (k < nyy-nb)) 
    {
        damper = damp1D[nb-(j-(nxx-nb))-1];
    }
    else if((i >= nb) && (i < nzz-nb) && (j >= nb) && (j < nxx-nb) && (k >= 0) && (k < nb)) 
    {
        damper = damp1D[k];
    }
    else if((i >= nb) && (i < nzz-nb) && (j >= nb) && (j < nxx-nb) && (k >= nyy-nb) && (k < nyy)) 
    {
        damper = damp1D[nb-(k-(nyy-nb))-1];
    }

    // 2D damping 
    else if((i >= nb) && (i < nzz-nb) && (j >= 0) && (j < nb) && (k >= 0) && (k < nb))
    {
        damper = damp2D[j + k*nb];
    }
    else if((i >= nb) && (i < nzz-nb) && (j >= nxx-nb) && (j < nxx) && (k >= 0) && (k < nb))
    {
        damper = damp2D[nb-(j-(nxx-nb))-1 + k*nb];
    }
    else if((i >= nb) && (i < nzz-nb) && (j >= 0) && (j < nb) && (k >= nyy-nb) && (k < nyy))
    {
        damper = damp2D[j + (nb-(k-(nyy-nb))-1)*nb];
    }
    else if((i >= nb) && (i < nzz-nb) && (j >= nxx-nb) && (j < nxx) && (k >= nyy-nb) && (k < nyy))
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

__global__ void compute_pressure(float * Txx, float * Tyy, float * Tzz, float * Pressure, int volsize)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index < volsize) Pressure[index] = (Txx[index] + Tyy[index] + Tzz[index]) / 3.0f;
}

__global__ void get_seismogram()
{



}















