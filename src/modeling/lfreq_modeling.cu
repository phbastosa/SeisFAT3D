# include "lfreq_modeling.cuh"

void lfreq_Modeling::display_progress()
{
    if (time_index % (nt / 10) == 0)
    {
        print_information();

        // stability and dispersion for 8E2T FDM operators

        std::cout<<"Time progress: " << floorf(100.0f * (float)(time_index+1) / (float)(nt)) <<" %\n\n";
    }
}

void lfreq_Modeling::define_cerjan_dampers()
{
    float * d1D = new float[nabc]();
    float * d2D = new float[nabc*nabc]();
    float * d3D = new float[nabc*nabc*nabc]();

    float factor = std::stof(catch_parameter("boundary_damper", file));

    for (int i = 0; i < nabc; i++) 
    {
        d1D[i] = expf(-powf(factor * (nabc - i), 2.0f));
    }

    for(int i = 0; i < nabc; i++) 
    {
        for (int j = 0; j < nabc; j++)
        {   
            d2D[j + i*nabc] += d1D[i]; // up to bottom
            d2D[i + j*nabc] += d1D[i]; // left to right
        }
    }

    for (int i  = 0; i < nabc; i++)
    {
        for(int j = 0; j < nabc; j++)
        {
            for(int k = 0; k < nabc; k++)
            {
                d3D[i + j*nabc + k*nabc*nabc] += d2D[i + j*nabc]; // XY plane
                d3D[i + j*nabc + k*nabc*nabc] += d2D[j + k*nabc]; // ZX plane
                d3D[i + j*nabc + k*nabc*nabc] += d2D[i + k*nabc]; // ZY plane
            }
        }
    }    

    for (int index = 0; index < nabc*nabc; index++)
        d2D[index] -= 1.0f;

    for (int index = 0; index < nabc*nabc*nabc; index++)
        d3D[index] -= 5.0f;    

	cudaMalloc((void**)&(damp1D), nabc*sizeof(float));
	cudaMalloc((void**)&(damp2D), nabc*nabc*sizeof(float));
	cudaMalloc((void**)&(damp3D), nabc*nabc*nabc*sizeof(float));

	cudaMemcpy(damp1D, d1D, nabc*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(damp2D, d2D, nabc*nabc*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(damp3D, d3D, nabc*nabc*nabc*sizeof(float), cudaMemcpyHostToDevice);

    delete[] d1D;
    delete[] d2D;
    delete[] d3D;    
}

void lfreq_Modeling::define_wavelet_signature()
{
    float * signal = new float[nt]();

    float pi = 4.0f*atanf(1.0f);

    float t0 = 2.0f*sqrtf(pi)/fmax;
    float fc = fmax/(3.0f * sqrtf(pi));

    for (int n = 0; n < nt; n++)
    {
        float td = n*dt - t0;

        float arg = pi*pi*pi*fmax*fmax*td*td;

        signal[n] = (1.0f - 2.0f*arg)*expf(-arg);
    }
    
    cudaMalloc((void**)&(wavelet), nt*sizeof(float));

    cudaMemcpy(wavelet, signal, nt*sizeof(float), cudaMemcpyHostToDevice);

    delete[] signal;
}

void lfreq_Modeling::define_grid_nodes_position()
{
    int * rx = new int[total_nodes]();
    int * ry = new int[total_nodes]();
    int * rz = new int[total_nodes]();

    for (int index = 0; index < total_nodes; index++)
    {
        rx[index] = (int)(geometry->nodes.x[index] / dx) + nbxl;
        ry[index] = (int)(geometry->nodes.y[index] / dy) + nbyl;
        rz[index] = (int)(geometry->nodes.z[index] / dz) + nbzu;
    }

	cudaMalloc((void**)&(grid_node_x), total_nodes*sizeof(int));
	cudaMalloc((void**)&(grid_node_y), total_nodes*sizeof(int));
	cudaMalloc((void**)&(grid_node_z), total_nodes*sizeof(int));

    cudaMemcpy(grid_node_x, rx, total_nodes*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(grid_node_y, ry, total_nodes*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(grid_node_z, rz, total_nodes*sizeof(int), cudaMemcpyHostToDevice);

    delete[] rx;
    delete[] ry;
    delete[] rz;
}

void lfreq_Modeling::set_outputs()
{   
    receiver_output_samples = nt*total_nodes;
    wavefield_output_samples = nPoints*total_snaps;

    if (export_receiver_output)
        receiver_output = new float[receiver_output_samples]();
    
    if (export_wavefield_output)
        wavefield_output = new float[wavefield_output_samples]();

    cudaMalloc((void**)&(seismogram), receiver_output_samples*sizeof(float));
}

void lfreq_Modeling::set_volumes()
{
    type_name = std::string("scalar");
    type_message = std::string("Constant density acoustic isotropic media");

    P = new float[volsize]();

    define_wavelet_signature();
    
    cudaMalloc((void**)&(d_Vp), volsize*sizeof(float));

    cudaMemcpy(d_Vp, V, volsize*sizeof(float), cudaMemcpyHostToDevice);

    cudaMalloc((void**)&(d_Upas), volsize*sizeof(float));
    cudaMalloc((void**)&(d_Upre), volsize*sizeof(float));
    cudaMalloc((void**)&(d_Ufut), volsize*sizeof(float));
}

void lfreq_Modeling::set_specifics()
{
    nt = std::stoi(catch_parameter("time_samples", file));
    dt = std::stof(catch_parameter("time_spacing", file));

    fmax = std::stof(catch_parameter("max_frequency", file));

    nabc = std::stoi(catch_parameter("boundary_samples", file));

    total_snaps = std::stoi(catch_parameter("total_snapshots", file));

    // free surface implementation

    nbxl = nabc; nbxr = nabc;
    nbyl = nabc; nbyr = nabc;
    nbzu = nabc; nbzd = nabc;    

    define_cerjan_dampers();
    define_grid_nodes_position();    
}

void lfreq_Modeling::initialization()
{
    cudaMemset(d_Upas, 0.0f, volsize*sizeof(float));
    cudaMemset(d_Upre, 0.0f, volsize*sizeof(float));
    cudaMemset(d_Ufut, 0.0f, volsize*sizeof(float));

    snap_index = 0;
}

void lfreq_Modeling::get_receiver_output()
{
    if (export_receiver_output)
    {
        cudaMemcpy(receiver_output, seismogram, nt*total_nodes*sizeof(float), cudaMemcpyDeviceToHost);

        receiver_output_file = receiver_output_folder + type_name + "_seismogram_Nsamples" + std::to_string(nt) + "_nRec" + std::to_string(total_nodes) + "_shot_" + std::to_string(shot_index+1) + ".bin";
    }
}

void lfreq_Modeling::get_wavefield_output()
{
    if (export_wavefield_output)
    {
        wavefield_output_file = wavefield_output_folder + type_name + "_snapshot_" + std::to_string(nz) + "x" + std::to_string(nx) + "x" + std::to_string(ny) + "_shot_" + std::to_string(shot_index+1) + "_Nsnaps" + std::to_string(total_snaps) + ".bin";
        
        if (snap_index < total_snaps)
        {
            if (time_index % (int)((float)(nt) / (float)(total_snaps)) == 0)
            {
                cudaMemcpy(P, d_Upre, volsize*sizeof(float), cudaMemcpyDeviceToHost);

                for (int index = 0; index < nPoints; index++)
                {
                    int y = (int) (index / (nx*nz));         
                    int x = (int) (index - y*nx*nz) / nz;    
                    int z = (int) (index - x*nz - y*nx*nz);  

                    wavefield_output[z + x*nz + y*nx*nz + snap_index*nPoints] = P[(z + nbzu) + (x + nbxl)*nzz + (y + nbyl)*nxx*nzz];
                }

                snap_index++;
            }
        }
    }
}

void lfreq_Modeling::get_seismogram()
{
    if (export_receiver_output)
    {
        int seismBlocks = (int)(total_nodes / threadsPerBlock) + 1;

        compute_seismogram<<<seismBlocks, threadsPerBlock>>>(seismogram, d_Upre, grid_node_x, grid_node_y, grid_node_z, total_nodes, nxx, nzz, nt, time_index);
    }
}

void lfreq_Modeling::forward_propagation()
{
    for (time_index = 0; time_index < nt; time_index++)
    {
        display_progress();

        compute_pressure<<<blocksPerGrid, threadsPerBlock>>>(d_Upre, d_Upas, d_Ufut, d_Vp, damp1D, damp2D, damp3D, wavelet, source_index, time_index, dx, dy, dz, dt, nxx, nyy, nzz, nabc);
        cudaDeviceSynchronize();

        update_pressure<<<blocksPerGrid, threadsPerBlock>>>(d_Upre, d_Upas, d_Ufut, volsize);
        cudaDeviceSynchronize();

        get_wavefield_output();
        get_seismogram();
    }   

    get_receiver_output();
}

void lfreq_Modeling::free_space()
{
    cudaFree(d_Vp);

    cudaFree(d_Upre);
    cudaFree(d_Upas);
    cudaFree(d_Ufut);
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

    float damper = get_boundary_damper(damp1D, damp2D, damp3D, i, j, k, nxx, nyy, nzz, nabc);
    
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

__global__ void compute_seismogram(float * seismogram, float * P, int * rx, int * ry, int * rz, int total_nodes, int nxx, int nzz, int nt, int time_id)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x; 

    if (index < total_nodes) 
        seismogram[time_id + index*nt] = P[rz[index] + rx[index]*nzz + ry[index]*nxx*nzz];
}

__device__ float get_boundary_damper(float * damp1D, float * damp2D, float * damp3D, int i, int j, int k, int nxx, int nyy, int nzz, int nabc)
{
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

    return damper;
}