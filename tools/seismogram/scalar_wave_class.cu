# include "scalar_wave_class.cuh"

void Acoustic::set_patameters()
{
    import_parameters();
    
    prepare_dampers();
    prepare_wavelet();
    prepare_geometry();
    prepare_velocity();
    prepare_wavefield();
}

void Acoustic::import_parameters()
{
    nx = std::stoi(catch_parameter("x_samples", file));
    ny = std::stoi(catch_parameter("y_samples", file));
    nz = std::stoi(catch_parameter("z_samples", file));
    nt = std::stoi(catch_parameter("time_samples", file));    

    dx = std::stof(catch_parameter("x_spacing", file));
    dy = std::stof(catch_parameter("y_spacing", file));
    dz = std::stof(catch_parameter("z_spacing", file));
    dt = std::stof(catch_parameter("time_spacing", file));
        
    nb = std::stoi(catch_parameter("boundary_samples", file)); 
}

void Acoustic::prepare_dampers()
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

void Acoustic::prepare_wavelet()
{
    float fmax = std::stof(catch_parameter("max_frequency", file));

    float pi = 4.0f * atanf(1.0f);
    float fc = fmax / (3.0f * sqrtf(pi)); 
    float to = 2.0f * sqrtf(pi) / fmax;

    float * ricker = new float[nt]();

    for (int n = 0; n < nt; n++)
    {        
        float arg = pi*powf(((n*dt - to)*fc*pi), 2.0f);
        
        ricker[n] = (1.0f - 2.0f*arg)*expf(-arg);    
    }

	cudaMalloc((void**)&(wavelet), nt*sizeof(float));
	cudaMemcpy(wavelet, ricker, nt*sizeof(float), cudaMemcpyHostToDevice);

    delete[] ricker;
}

void Acoustic::prepare_geometry()
{
    std::vector<std::string> shots;
    std::vector<std::string> nodes;
    std::vector<std::string> splitted;

    std::string shots_file = catch_parameter("shots_file", file);
    std::string nodes_file = catch_parameter("nodes_file", file);

    bool reciprocity = str2bool(catch_parameter("reciprocity", file));

    import_text_file(shots_file, shots);
    import_text_file(nodes_file, nodes);

    total_shots = shots.size();
    total_nodes = nodes.size();

    sx = new int[total_shots]();
    sy = new int[total_shots]();
    sz = new int[total_shots]();

    int * rx = new int[total_nodes]();
    int * ry = new int[total_nodes]();
    int * rz = new int[total_nodes]();

    for (int i = 0; i < total_shots; i++)
    {
        splitted = split(shots[i], ',');

        sx[i] = (int)(std::stof(splitted[0]) / dx) + nb;
        sy[i] = (int)(std::stof(splitted[1]) / dy) + nb;
        sz[i] = (int)(std::stof(splitted[2]) / dz) + nb;
    }    

    for (int i = 0; i < total_nodes; i++)
    {
        splitted = split(nodes[i], ',');

        rx[i] = (int)(std::stof(splitted[0]) / dx) + nb;
        ry[i] = (int)(std::stof(splitted[1]) / dy) + nb;
        rz[i] = (int)(std::stof(splitted[2]) / dz) + nb;
    }    

    std::vector<std::string>().swap(shots);
    std::vector<std::string>().swap(nodes);

    if (reciprocity)
    {
        std::swap(sx, rx);    
        std::swap(sy, ry);    
        std::swap(sz, rz);

        std::swap(total_shots, total_nodes);
    }

	cudaMalloc((void**)&(gx), total_nodes*sizeof(int));
	cudaMalloc((void**)&(gy), total_nodes*sizeof(int));
	cudaMalloc((void**)&(gz), total_nodes*sizeof(int));

	cudaMemcpy(gx, rx, total_nodes*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(gy, ry, total_nodes*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(gz, rz, total_nodes*sizeof(int), cudaMemcpyHostToDevice);

    delete[] rx;
    delete[] ry;
    delete[] rz;
}

void Acoustic::prepare_velocity()
{
    nxx = nx + 2*nb;
    nyy = ny + 2*nb;
    nzz = nz + 2*nb;

    model_size = nx*ny*nz;
    volume_size = nxx*nyy*nzz;

    std::string model_path = catch_parameter("model_path", file);

    vp = new float[model_size]();

    import_binary_float(model_path, vp, model_size);

    dtvp2 = new float[volume_size]();

    expand_model_boundary();

	cudaMalloc((void**)&(dtVp2), volume_size*sizeof(float));
    cudaMemcpy(dtVp2, dtvp2, volume_size*sizeof(float), cudaMemcpyHostToDevice);

    delete[] vp;
    delete[] dtvp2;
}

void Acoustic::expand_model_boundary()
{
    // Centering
    for (int z = nb; z < nzz - nb; z++)
    {
        for (int y = nb; y < nyy - nb; y++)
        {
            for (int x = nb; x < nxx - nb; x++)
            {
                dtvp2[z + x*nzz + y*nxx*nzz] = powf(dt*vp[(z - nb) + (x - nb)*nz + (y - nb)*nx*nz], 2.0f);
            }
        }
    }

    // Z direction
    for (int z = 0; z < nb; z++)
    {
        for (int y = nb; y < nyy - nb; y++)
        {
            for (int x = nb; x < nxx - nb; x++)
            {
                dtvp2[z + x*nzz + y*nxx*nzz] = powf(dt*vp[0 + (x - nb)*nz + (y - nb)*nx*nz], 2.0f);
                dtvp2[(nzz - z - 1) + x*nzz + y*nxx*nzz] = powf(dt*vp[(nz - 1) + (x - nb)*nz + (y - nb)*nx*nz], 2.0f);
            }
        }
    }

    // X direction
    for (int x = 0; x < nb; x++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int y = nb; y < nyy - nb; y++)
            {
                dtvp2[z + x*nzz + y*nxx*nzz] = dtvp2[z + nb*nzz + y*nxx*nzz];
                dtvp2[z + (nxx - x - 1)*nzz + y*nxx*nzz] = dtvp2[z + (nxx - nb - 1)*nzz + y*nxx*nzz];
            }
        }
    }

    // Y direction
    for (int y = 0; y < nb; y++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int x = 0; x < nxx; x++)
            {
                dtvp2[z + x*nzz + y*nxx*nzz] = dtvp2[z + x*nzz + nb*nxx*nzz];
                dtvp2[z + x*nzz + (nyy - y - 1)*nxx*nzz] = dtvp2[z + x*nzz + (nyy - nb - 1)*nxx*nzz];
            }
        }
    }
}

void Acoustic::prepare_wavefield()
{
    output_samples = nt * total_nodes;

    output_data = new float[output_samples]();

	cudaMalloc((void**)&(U), volume_size*sizeof(float));
	cudaMalloc((void**)&(Unew), volume_size*sizeof(float));
	cudaMalloc((void**)&(Uold), volume_size*sizeof(float));
    
	cudaMalloc((void**)&(seismogram), output_samples*sizeof(float));

    threadsPerBlock = 256;
    wavefieldBlocks = (int)(volume_size / threadsPerBlock);
    seismogramBlocks = (int)(output_samples / threadsPerBlock);
}

void Acoustic::info_message()
{
    if ((time_id - 1) % (nt / 10) == 0)
    {
        auto clear = system("clear");
            
        std::cout << "Model dimensions (z = " << (nz - 1)*dz << ", x = " << (nx - 1)*dx << ", y = " << (ny - 1)*dy << ") m\n\n";

        std::cout << "Shot "<< shot_id + 1 << " of " << total_shots;

        std::cout << " at position (z = " << (sz[shot_id] - nb)*dz << ", x = " 
                                          << (sx[shot_id] - nb)*dx << ", y = " 
                                          << (sy[shot_id] - nb)*dy << ") m\n\n";

        std::cout<<"Memory usage: \n";
        std::cout<<"RAM = "<<RAM<<" Mb\n";
        std::cout<<"GPU = "<<vRAM<<" Mb\n\n";

        std::cout<<"Modeling progress: " << 100.0f * (float)(time_id) / (float)(nt) << " % \n\n";
    } 
}

void Acoustic::initial_setup()
{
    get_RAM_usage();
    get_GPU_usage();

    cudaMemset(U, 0.0f, volume_size*sizeof(float));
    cudaMemset(Unew, 0.0f, volume_size*sizeof(float));
    cudaMemset(Uold, 0.0f, volume_size*sizeof(float));
    cudaMemset(seismogram, 0.0f, output_samples*sizeof(float));

    source_id = sz[shot_id] + sx[shot_id]*nzz + sy[shot_id]*nxx*nzz;    
}

void Acoustic::forward_solver()
{
    for (time_id = 0; time_id < nt; time_id++)
    {
        info_message();
        
        compute_wavefield<<<wavefieldBlocks, threadsPerBlock>>>(Unew, U, Uold, dtVp2, damp1D, damp2D, damp3D, wavelet, source_id, time_id, dx, dy, dz, nxx, nyy, nzz, nb);
        cudaDeviceSynchronize();
    
        compute_seismogram<<<seismogramBlocks, threadsPerBlock>>>(Unew, seismogram, gx, gy, gz, total_nodes, time_id, nt, nxx, nzz);
        cudaDeviceSynchronize();

        update_wavefield<<<wavefieldBlocks, threadsPerBlock>>>(Unew, U, Uold, volume_size);
        cudaDeviceSynchronize();
    }
}

void Acoustic::export_outputs()
{
    cudaMemcpy(output_data, seismogram, output_samples*sizeof(float), cudaMemcpyDeviceToHost);

    std::string data_path = "segy_data/synthetic_seismogram_" + std::to_string(nt) + "x" + std::to_string(total_nodes) + "_shot_" + std::to_string(shot_id+1) + ".bin";

    export_binary_float(data_path, output_data, output_samples);
}

void Acoustic::free_space()
{


}

void Acoustic::set_runtime()
{
    ti = std::chrono::system_clock::now();
}

void Acoustic::get_runtime()
{
    tf = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = tf - ti;

    std::ofstream runTimeFile("elapsedTime.txt",std::ios::in | std::ios::app);
    runTimeFile << "# Run Time [s]; RAM usage [MB]; GPU memory usage [MB]\n";
    runTimeFile << std::to_string(elapsed_seconds.count()) + "; " + std::to_string(RAM) + "; " + std::to_string(vRAM) + "\n";
    runTimeFile.close();

    std::cout<<"\nRun time: "<<elapsed_seconds.count()<<" s."<<std::endl;
}

void Acoustic::get_RAM_usage()
{
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    RAM = (int) (usage.ru_maxrss / 1024);
}

void Acoustic::get_GPU_usage()
{
	size_t freeMem, totalMem;
	cudaMemGetInfo(&freeMem, &totalMem);
    vRAM = (int) ((totalMem - freeMem) / (1024 * 1024));
    vRAM -= ivRAM;
}

void Acoustic::get_GPU_initMem()
{
	size_t freeMem, totalMem;
	cudaMemGetInfo(&freeMem, &totalMem);
    ivRAM = (int) ((totalMem - freeMem) / (1024 * 1024));
}

__global__ void update_wavefield(float * Unew, float * U, float * Uold, int volume_size)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (index < volume_size)
    {
        Uold[index] = U[index];        
        U[index] = Unew[index];
    }
}

__global__ void compute_wavefield(float * Unew, float * U, float * Uold, float * dtVp2, float * damp1D, float * damp2D, float * damp3D, float * wavelet, int source_id, int time_id, float dx, float dy, float dz, int nxx, int nyy, int nzz, int nb)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         // y direction
    int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
    int i = (int) (index - j*nzz - k*nxx*nzz); // z direction
    
    if (index == 0) U[source_id] += wavelet[time_id] / (dx*dy*dz);

    if((i >= 4) && (i < nzz-4) && (j >= 4) && (j < nxx-4) && (k >= 4) && (k < nyy-4)) 
    {
        float d2Px2 = (- 9.0f*(U[i + (j-4)*nzz + k*nxx*nzz] + U[i + (j+4)*nzz + k*nxx*nzz])
                   +   128.0f*(U[i + (j-3)*nzz + k*nxx*nzz] + U[i + (j+3)*nzz + k*nxx*nzz])
                   -  1008.0f*(U[i + (j-2)*nzz + k*nxx*nzz] + U[i + (j+2)*nzz + k*nxx*nzz])
                   +  8064.0f*(U[i + (j-1)*nzz + k*nxx*nzz] + U[i + (j+1)*nzz + k*nxx*nzz])
                   - 14350.0f*(U[i + j*nzz + k*nxx*nzz])) / (5040.0f*powf(dx, 2.0f));

        float d2Py2 = (- 9.0f*(U[i + j*nzz + (k-4)*nxx*nzz] + U[i + j*nzz + (k+4)*nxx*nzz])
                   +   128.0f*(U[i + j*nzz + (k-3)*nxx*nzz] + U[i + j*nzz + (k+3)*nxx*nzz])
                   -  1008.0f*(U[i + j*nzz + (k-2)*nxx*nzz] + U[i + j*nzz + (k+2)*nxx*nzz])
                   +  8064.0f*(U[i + j*nzz + (k-1)*nxx*nzz] + U[i + j*nzz + (k+1)*nxx*nzz])
                   - 14350.0f*(U[i + j*nzz + k*nxx*nzz])) / (5040.0f*powf(dy,2.0f));

        float d2Pz2 = (- 9.0f*(U[(i-4) + j*nzz + k*nxx*nzz] + U[(i+4) + j*nzz + k*nxx*nzz])
                   +   128.0f*(U[(i-3) + j*nzz + k*nxx*nzz] + U[(i+3) + j*nzz + k*nxx*nzz])
                   -  1008.0f*(U[(i-2) + j*nzz + k*nxx*nzz] + U[(i+2) + j*nzz + k*nxx*nzz])
                   +  8064.0f*(U[(i-1) + j*nzz + k*nxx*nzz] + U[(i+1) + j*nzz + k*nxx*nzz])
                   - 14350.0f*(U[i + j*nzz + k*nxx*nzz])) / (5040.0f*powf(dz,2.0f));
        
        Unew[index] = dtVp2[index] * (d2Px2 + d2Py2 + d2Pz2) + 2.0f*U[index] - Uold[index]; 
    }

    float damper = get_boundary_damper(damp1D, damp2D, damp3D, i, j, k, nxx, nyy, nzz, nb);
    
    if (index < nxx*nyy*nzz)
    {
        U[index] *= damper;
        Uold[index] *= damper;    
        Unew[index] *= damper;
    }
}

__global__ void compute_seismogram(float * Unew, float * seismogram, int * gx, int * gy, int * gz, int total_nodes, int time_id, int nt, int nxx, int nzz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x; 

    if (index < total_nodes) 
        seismogram[time_id + index*nt] = Unew[gz[index] + gx[index]*nzz + gy[index]*nxx*nzz];
}

__device__ float get_boundary_damper(float * damp1D, float * damp2D, float * damp3D, int i, int j, int k, int nxx, int nyy, int nzz, int nb)
{
    if ((i >= nb) && (i < nzz-nb) && (j >= nb) && (j < nxx-nb) && (k >= nb) && (k < nyy-nb))  
        return 1.0f;

    // 1D damping
    
    if((i < nb) && (j >= nb) && (j < nxx-nb) && (k >= nb) && (k < nyy-nb)) 
        return damp1D[i];
    
    if((i >= nzz-nb) && (i < nzz) && (j >= nb) && (j < nxx-nb) && (k >= nb) && (k < nyy-nb)) 
        return damp1D[nb-(i-(nzz-nb))-1];
    
    if((i >= nb) && (i < nzz-nb) && (j >= 0) && (j < nb) && (k >= nb) && (k < nyy-nb)) 
        return damp1D[j];
    
    if((i >= nb) && (i < nzz-nb) && (j >= nxx-nb) && (j < nxx) && (k >= nb) && (k < nyy-nb)) 
        return damp1D[nb-(j-(nxx-nb))-1];
    
    if((i >= nb) && (i < nzz-nb) && (j >= nb) && (j < nxx-nb) && (k >= 0) && (k < nb)) 
        return damp1D[k];
    
    if((i >= nb) && (i < nzz-nb) && (j >= nb) && (j < nxx-nb) && (k >= nyy-nb) && (k < nyy)) 
        return damp1D[nb-(k-(nyy-nb))-1];
    
    // 2D damping 
    if((i >= nb) && (i < nzz-nb) && (j >= 0) && (j < nb) && (k >= 0) && (k < nb))
        return damp2D[j + k*nb];
    
    if((i >= nb) && (i < nzz-nb) && (j >= nxx-nb) && (j < nxx) && (k >= 0) && (k < nb))
        return damp2D[nb-(j-(nxx-nb))-1 + k*nb];
    
    if((i >= nb) && (i < nzz-nb) && (j >= 0) && (j < nb) && (k >= nyy-nb) && (k < nyy))
        return damp2D[j + (nb-(k-(nyy-nb))-1)*nb];
    
    if((i >= nb) && (i < nzz-nb) && (j >= nxx-nb) && (j < nxx) && (k >= nyy-nb) && (k < nyy))
        return damp2D[nb-(j-(nxx-nb))-1 + (nb-(k-(nyy-nb))-1)*nb];
    
    if((i >= 0) && (i < nb) && (j >= nb) && (j < nxx-nb) && (k >= 0) && (k < nb))
        return damp2D[i + k*nb];
    
    if((i >= nzz-nb) && (i < nzz) && (j >= nb) && (j < nxx-nb) && (k >= 0) && (k < nb))
        return damp2D[nb-(i-(nzz-nb))-1 + k*nb];
    
    if((i >= 0) && (i < nb) && (j >= nb) && (j < nxx-nb) && (k >= nyy-nb) && (k < nyy))
        return damp2D[i + (nb-(k-(nyy-nb))-1)*nb];
    
    if((i >= nzz-nb) && (i < nzz) && (j >= nb) && (j < nxx-nb) && (k >= nyy-nb) && (k < nyy))
        return damp2D[nb-(i-(nzz-nb))-1 + (nb-(k-(nyy-nb))-1)*nb];
    
    if((i >= 0) && (i < nb) && (j >= 0) && (j < nb) && (k >= nb) && (k < nyy-nb))
        return damp2D[i + j*nb];
    
    if((i >= nzz-nb) && (i < nzz) && (j >= 0) && (j < nb) && (k >= nb) && (k < nyy-nb))
        return damp2D[nb-(i-(nzz-nb))-1 + j*nb];
    
    if((i >= 0) && (i < nb) && (j >= nxx-nb) && (j < nxx) && (k >= nb) && (k < nyy-nb))
        return damp2D[i + (nb-(j-(nxx-nb))-1)*nb];
    
    if((i >= nzz-nb) && (i < nzz) && (j >= nxx-nb) && (j < nxx) && (k >= nb) && (k < nyy-nb))
        return damp2D[nb-(i-(nzz-nb))-1 + (nb-(j-(nxx-nb))-1)*nb];
    
    // 3D damping
    if((i >= 0) && (i < nb) && (j >= 0) && (j < nb) && (k >= 0) && (k < nb))
        return damp3D[i + j*nb + k*nb*nb];
    
    if((i >= nzz-nb) && (i < nzz) && (j >= 0) && (j < nb) && (k >= 0) && (k < nb))
        return damp3D[nb-(i-(nzz-nb))-1 + j*nb + k*nb*nb];
    
    if((i >= 0) && (i < nb) && (j >= nxx-nb) && (j < nxx) && (k >= 0) && (k < nb))
        return damp3D[i + (nb-(j-(nxx-nb))-1)*nb + k*nb*nb];
    
    if((i >= 0) && (i < nb) && (j >= 0) && (j < nb) && (k >= nyy-nb) && (k < nyy))
        return damp3D[i + j*nb + (nb-(k-(nyy-nb))-1)*nb*nb];
    
    if((i >= nzz-nb) && (i < nzz) && (j >= nxx-nb) && (j < nxx) && (k >= 0) && (k < nb))
        return damp3D[nb-(i-(nzz-nb))-1 + (nb-(j-(nxx-nb))-1)*nb + k*nb*nb];
    
    if((i >= nzz-nb) && (i < nzz) && (j >= 0) && (j < nb) && (k >= nyy-nb) && (k < nyy))
        return damp3D[nb-(i-(nzz-nb))-1 + j*nb + (nb-(k-(nyy-nb))-1)*nb*nb];
    
    if((i >= 0) && (i < nb) && (j >= nxx-nb) && (j < nxx) && (k >= nyy-nb) && (k < nyy))
        return damp3D[i + (nb-(j-(nxx-nb))-1)*nb + (nb-(k-(nyy-nb))-1)*nb*nb];
    
    if((i >= nzz-nb) && (i < nzz) && (j >= nxx-nb) && (j < nxx) && (k >= nyy-nb) && (k < nyy))
        return damp3D[nb-(i-(nzz-nb))-1 + (nb-(j-(nxx-nb))-1)*nb + (nb-(k-(nyy-nb))-1)*nb*nb];
}
