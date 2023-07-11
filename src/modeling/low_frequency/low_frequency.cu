# include "low_frequency.cuh"

void Low_Frequency::set_specifications()
{
    pi = 4.0f * atanf(1.0f);  

    tmax = std::stoi(catch_parameter("total_time", file));
    dt = std::stof(catch_parameter("time_spacing", file));
    nt = (int)(ceilf(tmax / dt)) + 1; 

    nsnap = std::stoi(catch_parameter("n_snapshots", file));
    dsnap = (int)(nt / (nsnap-1));

    free_surface = str2bool(catch_parameter("free_surface", file));

    fmax = std::stof(catch_parameter("max_frequency", file));
    tlag = std::stof(catch_parameter("wavelet_shift", file));
    amp = std::stof(catch_parameter("max_amplitude", file));
    fc = fmax / (3.0f * sqrtf(pi));

    set_model_boundaries();
    set_model_parameters();
    set_gridded_geometry();
    set_wavefields();

    set_wavelet();
    set_dampers();
    set_outputs();
}

void Low_Frequency::set_model_boundaries()
{
    nbxl = nb; nbxr = nb;
    nbyl = nb; nbyr = nb;
    nbzu = nb; nbzd = nb;

    if (free_surface) nbzu = 4; // FDM operators: 8E2T

    nxx = nx + nbxl + nbxr;
    nyy = ny + nbyl + nbyr;
    nzz = nz + nbzu + nbzd;

    volsize = nxx*nyy*nzz;

    blocksPerGrid = (int)(volsize / threadsPerBlock);
}

void Low_Frequency::set_gridded_geometry()
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

void Low_Frequency::set_dampers()
{
    float * d1D = new float[nb]();
    float * d2D = new float[nb*nb]();
    float * d3D = new float[nb*nb*nb]();

    float factor = std::stof(catch_parameter("boundary_damper", file));

    for (int i = 0; i < nb; i++) 
    {
        d1D[i] = expf(-powf(factor * (nb - i), 2.0f));
    }

    for(int i = 0; i < nb; i++) 
    {
        for (int j = 0; j < nb; j++)
        {   
            d2D[j + i*nb] += d1D[i]; // up to bottom
            d2D[i + j*nb] += d1D[i]; // left to right
        }
    }

    for (int i  = 0; i < nb; i++)
    {
        for(int j = 0; j < nb; j++)
        {
            for(int k = 0; k < nb; k++)
            {
                d3D[i + j*nb + k*nb*nb] += d2D[i + j*nb]; // XY plane
                d3D[i + j*nb + k*nb*nb] += d2D[j + k*nb]; // ZX plane
                d3D[i + j*nb + k*nb*nb] += d2D[i + k*nb]; // ZY plane
            }
        }
    }    

    for (int index = 0; index < nb*nb; index++)
        d2D[index] -= 1.0f;

    for (int index = 0; index < nb*nb*nb; index++)
        d3D[index] -= 5.0f;    

	cudaMalloc((void**)&(damp1D), nb*sizeof(float));
	cudaMalloc((void**)&(damp2D), nb*nb*sizeof(float));
	cudaMalloc((void**)&(damp3D), nb*nb*nb*sizeof(float));

	cudaMemcpy(damp1D, d1D, nb*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(damp2D, d2D, nb*nb*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(damp3D, d3D, nb*nb*nb*sizeof(float), cudaMemcpyHostToDevice);

    delete[] d1D;
    delete[] d2D;
    delete[] d3D;
}

void Low_Frequency::set_outputs()
{
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

void Low_Frequency::build_outputs()
{
    if (export_receiver_output)
    {
        cudaMemcpy(receiver_output, seismogram, receiver_output_samples*sizeof(float), cudaMemcpyDeviceToHost);
        receiver_output_file = receiver_output_folder + modeling_method + "_pressure_" + std::to_string(nt) + "x" + std::to_string(geometry->nodes.total) + "_shot_" + std::to_string(shot_id+1) + ".bin";
    }

    if (export_wavefield_output)
        wavefield_output_file = wavefield_output_folder + modeling_method + "_pressure_" + std::to_string(nz) + "x" + std::to_string(nx) + "x" + std::to_string(ny) + "_" + std::to_string(nsnap) + "_shot_" + std::to_string(shot_id+1) + ".bin";
}

void Low_Frequency::show_progress()
{
    if (time_id % (nt / 10) == 0) 
        info_message();
}

void Low_Frequency::get_snapshots()
{
    if (export_wavefield_output)
    {
        if (time_id % dsnap == 0)
        {
            cudaMemcpy(snapshot, Pressure, volsize*sizeof(float), cudaMemcpyDeviceToHost);       

            # pragma omp parallel for 
            for (int index = 0; index < nPoints; index++)
            {
                int y = (int) (index / (nx*nz));         
                int x = (int) (index - y*nx*nz) / nz;    
                int z = (int) (index - x*nz - y*nx*nz);  

                wavefield_output[z + x*nz + y*nx*nz + isnap*nPoints] = snapshot[(z + nbzu) + (x + nbxl)*nzz + (y + nbyl)*nxx*nzz];
            }

            isnap++;
        }
    }
}

void Low_Frequency::get_seismogram()
{
    if (export_receiver_output)
    {
        int seismBlocks = (int)(total_nodes / threadsPerBlock) + 1;

        compute_seismogram<<<seismBlocks,threadsPerBlock>>>(Pressure,seismogram,grid_node_x,grid_node_y,grid_node_z,total_nodes,time_id,nt,nxx,nzz);
    }    
}

__global__ void compute_seismogram(float * Pressure, float * seismogram, int * rx, int * ry, int * rz, int total_nodes, int time_id, int nt, int nxx, int nzz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x; 

    if (index < total_nodes) 
        seismogram[time_id + index*nt] = Pressure[rz[index] + rx[index]*nzz + ry[index]*nxx*nzz];
}

__device__ float get_boundary_damper(int i, int j, int k, int nxx, int nyy, int nzz, int nb, int nbzu)
{
    float damper = 1.0f;

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

    return damper;
}