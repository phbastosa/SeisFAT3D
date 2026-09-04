# include "migration.cuh"

void Migration::set_parameters()
{
    nt = std::stoi(catch_parameter("time_samples", parameters));
    dt = std::stof(catch_parameter("time_spacing", parameters));
    da = std::stof(catch_parameter("mig_angle_spacing", parameters));

    nsx = std::stoi(catch_parameter("nsx", parameters));
    nsy = std::stoi(catch_parameter("nsy", parameters));
    nrx = std::stoi(catch_parameter("nrx", parameters));
    nry = std::stoi(catch_parameter("nry", parameters));

    max_angle = std::stof(catch_parameter("mig_max_angle", parameters));
    max_it = std::stoi(catch_parameter("mig_max_iteration", parameters));
    residuo_folder = catch_parameter("mig_residuo_folder", parameters);

    fmax = std::stof(catch_parameter("max_frequency", parameters));

    aperture = std::stof(catch_parameter("mig_aperture", parameters));

    max_offset = std::stof(catch_parameter("mig_max_offset", parameters));

    input_data_folder = catch_parameter("mig_data_folder", parameters);
    input_data_prefix = catch_parameter("mig_data_prefix", parameters);
    
    tables_folder = catch_parameter("mig_tables_folder", parameters);
    seismic_folder = catch_parameter("mig_images_folder", parameters);

    anisotropy = str2bool(catch_parameter("anisotropy", parameters));

    if (anisotropy) modeling = new Eikonal_ANI(); 
               else modeling = new Eikonal_ISO();

    modType = anisotropy ? "anisotropic" : "isotropic";

    modeling->parameters = parameters;
    
    modeling->set_geometry();
    
    set_interpolation();
    set_anisotropy();
    set_slowness();
    set_wavelet();

    modeling->set_eikonal();
    modeling->set_conditions();
    
    set_migration();

    volBytes = modeling->volsize*sizeof(float);

    h_Ts = new float[nsy*modeling->volsize]();
    h_Tr = new float[nrx*modeling->volsize]();

    h_xsrc = new float[nsy]();
    h_ysrc = new float[nsy]();

    h_xrec = new float[nrx]();
    h_yrec = new float[nrx]();

    h_data = new float[nt*nrx*nsy]();
    h_model = new float[m_samples]();

    seismic = new float[nsy*nt*modeling->geometry->nrec]();

    cudaMalloc((void**)&d_Ts, nsy*volBytes);
    cudaMalloc((void**)&d_Tr, nrx*volBytes);

    cudaMalloc((void**)&d_xsrc, nsy*sizeof(float));
    cudaMalloc((void**)&d_ysrc, nsy*sizeof(float));

    cudaMalloc((void**)&d_xrec, nrx*sizeof(float));
    cudaMalloc((void**)&d_yrec, nrx*sizeof(float));

    cudaMalloc((void**)&d_data, nt*nrx*nsy*sizeof(float));
    cudaMalloc((void**)&d_model, m_samples*sizeof(float));
}

void Migration::set_interpolation()
{
    old_nx = std::stoi(catch_parameter("x_samples", parameters));
    old_ny = std::stoi(catch_parameter("y_samples", parameters));
    old_nz = std::stoi(catch_parameter("z_samples", parameters));

    old_dh = std::stof(catch_parameter("model_spacing", parameters));
    new_dh = std::stof(catch_parameter("cubic_spacing", parameters)); 

    new_nx = (int)((old_nx-1)*old_dh/new_dh)+1;
    new_ny = (int)((old_ny-1)*old_dh/new_dh)+1;
    new_nz = (int)((old_nz-1)*old_dh/new_dh)+1;

    old_nPoints = old_nx*old_ny*old_nz;
    new_nPoints = new_nx*new_ny*new_nz;

    modeling->nb = 3;

    modeling->dh = new_dh;

    modeling->nx = new_nx;
    modeling->ny = new_ny;
    modeling->nz = new_nz;

    modeling->nxx = modeling->nx + 2*modeling->nb;
    modeling->nyy = modeling->ny + 2*modeling->nb;
    modeling->nzz = modeling->nz + 2*modeling->nb;

    modeling->nPoints = new_nPoints;    
    modeling->volsize = modeling->nxx*modeling->nyy*modeling->nzz;

    modeling->nBlocks = (int)((modeling->volsize + NTHREADS - 1) / NTHREADS);

    nBlocks = (int)((old_nPoints + NTHREADS - 1) / NTHREADS);
}

void Migration::perform_cubic(float * input, float * output)
{
    float P[4][4][4];

    if ((new_nx == old_nx) && (new_ny == old_ny) && (new_nz == old_nz))
    {
        std::swap(output, input);
    }
    else
    {
        for (int index = 0; index < new_nPoints; index++)
        {
            int new_k = (int) (index / (new_nx*new_nz));         
            int new_j = (int) (index - new_k*new_nx*new_nz) / new_nz;    
            int new_i = (int) (index - new_j*new_nz - new_k*new_nx*new_nz); 
            
            float z = (float)(new_i) * new_dh;
            float x = (float)(new_j) * new_dh;
            float y = (float)(new_k) * new_dh;
            
            float z0 = floorf(z / old_dh) * old_dh;
            float x0 = floorf(x / old_dh) * old_dh;
            float y0 = floorf(y / old_dh) * old_dh;
            
            float zd = (z - z0) / old_dh;
            float xd = (x - x0) / old_dh;
            float yd = (y - y0) / old_dh;

            int old_i = (int)(z / old_dh); 
            int old_j = (int)(x / old_dh);   
            int old_k = (int)(y / old_dh);         

            if ((new_i > 0) && (new_i < new_nz-1) && (new_j > 0) && (new_j < new_nx-1) && (new_k > 0) && (new_k < new_ny-1))
            {
                for (int pIdx = 0; pIdx < 4; pIdx++)
                    for (int pIdy = 0; pIdy < 4; pIdy++)
                        for (int pIdz = 0; pIdz < 4; pIdz++)
                            P[pIdx][pIdy][pIdz] = input[(old_i + pIdz - 1) + (old_j + pIdx - 1)*old_nz + (old_k + pIdy - 1)*old_nx*old_nz];

                output[new_i + new_j*new_nz + new_k*new_nx*new_nz] = cubic3d(P, xd, yd, zd);
            }        
        }

        for (int i = 0; i < new_nz; i++)
        {
            for (int j = 0; j < new_nx; j++)
            {
                int beg = 0;
                int end = new_ny-1;

                output[i + j*new_nz + beg*new_nx*new_nz] = output[i + j*new_nz + (beg+1)*new_nx*new_nz];
                output[i + j*new_nz + end*new_nx*new_nz] = output[i + j*new_nz + (end-1)*new_nx*new_nz];
            }
        
            for (int k = 0; k < new_ny; k++)
            {
                int beg = 0;
                int end = new_nx-1;

                output[i + beg*new_nz + k*new_nx*new_nz] = output[i + (beg+1)*new_nz + k*new_nx*new_nz];
                output[i + end*new_nz + k*new_nx*new_nz] = output[i + (end-1)*new_nz + k*new_nx*new_nz];
            }    
        }    

        for (int j = 0; j < new_nx; j++)
        {
            for (int k = 0; k < new_ny; k++)
            {
                int beg = 0;
                int end = new_nz-1;

                output[beg + j*new_nz + k*new_nx*new_nz] = output[(beg+1) + j*new_nz + k*new_nx*new_nz];
                output[end + j*new_nz + k*new_nx*new_nz] = output[(end-1) + j*new_nz + k*new_nx*new_nz];
            }
        }
    }
}

void Migration::set_anisotropy()
{
    if (anisotropy)
    {
        modeling->Cijkl_folder = catch_parameter("Cijkl_folder", parameters);

        float * Cin = new float[old_nPoints]();
        float * Cout = new float[new_nPoints]();

        int n = 6;

        for (int i = 0; i < n; i++)
        {   
            for (int j = i; j < n; j++)
            {
                std::string Cij = "C" + std::to_string(i+1) + std::to_string(j+1);
            
                import_binary_float(modeling->Cijkl_folder + Cij + ".bin", Cin, old_nPoints);
                perform_cubic(Cin, Cout);
                export_binary_float(modeling->Cijkl_folder + "cubic_" + Cij + ".bin", Cout, new_nPoints); 
            }
        }

        modeling->Cijkl_folder = modeling->Cijkl_folder + "cubic_";
   
        delete[] Cin;
        delete[] Cout;
    }
}

void Migration::set_slowness()
{
    float * old_s = new float[old_nPoints]();
    float * new_s = new float[new_nPoints]();

    std::string slowness_file = catch_parameter("slowness_file", parameters);

    import_binary_float(slowness_file, old_s, old_nPoints);

    perform_cubic(old_s, new_s);

    modeling->S = new float[modeling->volsize]();

    modeling->expand_boundary(new_s, modeling->S);

    delete[] old_s;
    delete[] new_s;
}

void Migration::set_wavelet()
{
    float t0 = 2.0f*sqrtf(M_PI) / fmax;
    float fc = fmax / (3.0f * sqrtf(M_PI));

    nw = 2*((int)((t0 - 0.5f*dt) / dt));

    wavelet = new float[nw]();

    for (int wId = 0; wId < nw; wId++)
    {
        float td = wId*dt - t0;

        float arg = M_PI*M_PI*M_PI*fc*fc*td*td;

        wavelet[wId] = (1.0f - 2.0f*arg)*expf(-arg);
    }
}

void Migration::set_CMP_gathers()
{
    nang = (int)(max_angle / da) + 1;

    ds = modeling->geometry->ysrc[1] - modeling->geometry->ysrc[0];
    dr = modeling->geometry->yrec[1] - modeling->geometry->yrec[0];

    dCMP = 0.5f*min(ds,dr);

    minCMPx = 1e6f; minCMPy = 1e6f;
    maxCMPy =-1e6f; maxCMPx =-1e6f;  

    for (int srcId = 0; srcId < modeling->geometry->nsrc; srcId++)
    {
        for (int recId = 0; recId < modeling->geometry->nrec; recId++)
        {
            float CMPx = 0.5f*(modeling->geometry->xsrc[srcId] + 
                               modeling->geometry->xrec[recId]);

            float CMPy = 0.5f*(modeling->geometry->ysrc[srcId] + 
                               modeling->geometry->yrec[recId]);

            minCMPx = (CMPx < minCMPx) ? CMPx : minCMPx;
            minCMPy = (CMPy < minCMPy) ? CMPy : minCMPy;

            maxCMPx = (CMPx > maxCMPx) ? CMPx : maxCMPx;
            maxCMPy = (CMPy > maxCMPy) ? CMPy : maxCMPy;
        }
    }

    nCMPx = (int)((maxCMPx - minCMPx) / dCMP) + 1;
    nCMPy = (int)((maxCMPy - minCMPy) / dCMP) + 1;
    
    nCMP = nCMPx * nCMPy;
}

void Migration::prepare_convolution()
{
    nfft = nextpow2(nt + nw - 1);
    nfreq = nfft/2 + 1;

    time_trace = (double *) fftw_malloc(nfft*sizeof(double));
    time_wavelet = (double *) fftw_malloc(nfft*sizeof(double));

    freq_trace = (fftw_complex *) fftw_malloc(nfreq*sizeof(fftw_complex));
    freq_wavelet = (fftw_complex *) fftw_malloc(nfreq*sizeof(fftw_complex));

    trace_forward_plan = fftw_plan_dft_r2c_1d(nfft, time_trace, freq_trace, FFTW_ESTIMATE);
    trace_inverse_plan = fftw_plan_dft_c2r_1d(nfft, freq_trace, time_trace, FFTW_ESTIMATE);
    wavelet_forward_plan = fftw_plan_dft_r2c_1d(nfft, time_wavelet, freq_wavelet, FFTW_ESTIMATE);

    for (int tId = 0; tId < nfft; tId++)
    {
        time_trace[tId] = 0.0;
        time_wavelet[tId] = (tId < nw) ? (double)wavelet[tId] : 0.0;
    }

    fftw_execute(wavelet_forward_plan);
}

void Migration::get_src_travel_times()
{
    modeling->keyword = "source";
    modeling->total = std::to_string(modeling->geometry->nsrc); 

    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nsrc; modeling->srcId++)
    {
        modeling->sx = modeling->geometry->xsrc[modeling->srcId];
        modeling->sy = modeling->geometry->ysrc[modeling->srcId];
        modeling->sz = modeling->geometry->zsrc[modeling->srcId];

        modeling->current = std::to_string(modeling->srcId+1);
        
        modeling->xpos = format1Decimal(modeling->sx);
        modeling->ypos = format1Decimal(modeling->sy);
        modeling->zpos = format1Decimal(modeling->sz);

        modeling->current_operation = "Computing " + modeling->keyword + " travel time matrices in " + modType + " media";

        modeling->show_information();
        modeling->time_propagation();
        
        cudaMemcpy(modeling->T, modeling->d_T, modeling->volsize*sizeof(float), cudaMemcpyDeviceToHost);

        export_binary_float(tables_folder + "eikonal_src_" + std::to_string(modeling->srcId+1) + ".bin", modeling->T, modeling->volsize);
    }
}

void Migration::get_rec_travel_times()
{
    modeling->keyword = "receiver";
    modeling->total = std::to_string(modeling->geometry->nrec); 
    
    for (modeling->recId = 0; modeling->recId < modeling->geometry->nrec; modeling->recId++)
    {
        modeling->sx = modeling->geometry->xrec[modeling->recId];
        modeling->sy = modeling->geometry->yrec[modeling->recId];
        modeling->sz = modeling->geometry->zrec[modeling->recId];

        modeling->current = std::to_string(modeling->recId+1);
        
        modeling->xpos = format1Decimal(modeling->sx);
        modeling->ypos = format1Decimal(modeling->sy);
        modeling->zpos = format1Decimal(modeling->sz);

        modeling->current_operation = "Computing " + modeling->keyword + " travel time matrices in " + modType + " media";

        modeling->show_information();
        modeling->time_propagation();

        cudaMemcpy(modeling->T, modeling->d_T, modeling->volsize*sizeof(float), cudaMemcpyDeviceToHost);

        export_binary_float(tables_folder + "eikonal_rec_" + std::to_string(modeling->recId+1) + ".bin", modeling->T, modeling->volsize);
    }
}

void Migration::set_cross_spread_data()
{
    for (int src_csIdy = 0; src_csIdy < nsy; src_csIdy++)
    {
        int d_src_offset = src_csIdy*nt*nrx;
        int h_src_offset = src_csIdy*nt*modeling->geometry->nrec;

        for (int rec_csIdx = 0; rec_csIdx < nrx; rec_csIdx++)
        {            
            int d_rec_offset = rec_csIdx*nt;
            int h_rec_offset = (rec_csIdy + rec_csIdx*nry)*nt;

            int hst_base = h_rec_offset + h_src_offset; 
            int dvc_base = d_rec_offset + d_src_offset;

            for (int tId = 0; tId < nt; tId++)
                h_data[tId + dvc_base] = seismic[tId + hst_base];
        }
    }        
}

void Migration::set_src_line_components()
{
    for (int src_csIdy = 0; src_csIdy < nsy; src_csIdy++)
    {
        int srcId = src_csIdy + src_csIdx * nsy;

        float * target_seismic = seismic + src_csIdy * nt * modeling->geometry->nrec;
        float * target_eikonal = h_Ts + src_csIdy * modeling->volsize;

        std::string seismic_path = input_data_folder + input_data_prefix + std::to_string(srcId + 1) + ".bin";
        std::string eikonal_path = tables_folder + "eikonal_src_" + std::to_string(srcId + 1) + ".bin";

        import_binary_float(seismic_path, target_seismic, nt*modeling->geometry->nrec);
        import_binary_float(eikonal_path, target_eikonal, modeling->volsize);

        h_xsrc[src_csIdy] = modeling->geometry->xsrc[srcId];
        h_ysrc[src_csIdy] = modeling->geometry->ysrc[srcId];
    }

    modeling->total   = std::to_string(nsx);
    modeling->current = std::to_string(src_csIdx+1);

    modeling->zpos = format1Decimal(modeling->geometry->zsrc[src_csIdx*nsy]);
    modeling->xpos = format1Decimal(modeling->geometry->xsrc[src_csIdx*nsy]);
    modeling->ypos = format1Decimal(modeling->geometry->ysrc[    0   + src_csIdx*nsy]) + " - " +
                     format1Decimal(modeling->geometry->ysrc[(nsy-1) + src_csIdx*nsy]);
}

void Migration::set_rec_line_components()
{
    for (int rec_csIdx = 0; rec_csIdx < nrx; rec_csIdx++) 
    {
        int recId = rec_csIdy + rec_csIdx*nry;

        float * target = h_Tr + rec_csIdx*modeling->volsize;
        
        std::string path = tables_folder + "eikonal_rec_" + std::to_string(recId+1) + ".bin";
        
        import_binary_float(path, target, modeling->volsize);

        h_xrec[rec_csIdx] = modeling->geometry->xrec[recId];
        h_yrec[rec_csIdx] = modeling->geometry->yrec[recId];
    }

    zpos = format1Decimal(modeling->geometry->zrec[rec_csIdy]);
    ypos = format1Decimal(modeling->geometry->yrec[rec_csIdy]);
    xpos = format1Decimal(modeling->geometry->xrec[rec_csIdy + 0*nry]) + " - " +
            format1Decimal(modeling->geometry->xrec[rec_csIdy + (nrx-1)*nry]);    
}

void Migration::adjoint_convolution()
{
    for (int csId = 0; csId < nrx * nsy; csId++)
    {
        for (int tId = 0; tId < nfft; tId++)
        {
            if (tId < nt)
            {
                time_trace[tId] = (double)h_data[tId + csId*nt];
                h_data[tId + csId*nt] = 0.0f;
            }
            else 
                time_trace[tId] = 0.0;    
        }
        
        fftw_execute(trace_forward_plan);

        for (int fId = 0; fId < nfreq; fId++)
        {
            double a_re = freq_trace[fId][0];
            double a_im = freq_trace[fId][1];
            double b_re = freq_wavelet[fId][0];
            double b_im = freq_wavelet[fId][1];

            freq_trace[fId][0] = a_re*b_re + a_im*b_im;  
            freq_trace[fId][1] = a_im*b_re - a_re*b_im;  
        }

        fftw_execute(trace_inverse_plan);

        for (int tId = nw/2 + nw/5; tId < nt; tId++)
            h_data[tId + csId*nt] = (float)(time_trace[tId - nw/2 - nw/5]);
    }
}

void Migration::forward_convolution()
{
    for (int csId = 0; csId < nrx * nsy; csId++)
    {
        for (int tId = 0; tId < nfft; tId++)
        {
            if (tId < nt)
            {
                time_trace[tId] = (double)h_data[tId + csId*nt];
                h_data[tId + csId*nt] = 0.0f;
            }
            else 
                time_trace[tId] = 0.0;    
        }

        fftw_execute(trace_forward_plan);

        for (int fId = 0; fId < nfreq; fId++)
        {
            double a_re = freq_trace[fId][0];
            double a_im = freq_trace[fId][1];

            double b_re = freq_wavelet[fId][0];
            double b_im = freq_wavelet[fId][1];

            freq_trace[fId][0] = a_re * b_re - a_im * b_im;
            freq_trace[fId][1] = a_re * b_im + a_im * b_re;
        }

        fftw_execute(trace_inverse_plan);

        for (int tId = 0; tId < nt; tId++)
            h_data[tId + csId*nt] = (float)(time_trace[tId + nw/2]);
    }        
}

void Migration::show_information()
{
    auto clear = system("clear");
    
    std::cout << "-------------------------------------------------------------------\n";
    std::cout << " \033[34mSeisFAT3D\033[0;0m --------------------------------------------------------\n";
    std::cout << "-------------------------------------------------------------------\n\n";

    std::cout << "Model dimensions: (z = " << (old_nz - 1)*old_dh << 
                                  ", x = " << (old_nx - 1)*old_dh <<
                                  ", y = " << (old_ny - 1)*old_dh << ") m\n\n";

    std::cout << "Running Cross Spread for xline " << modeling->current << " of " << modeling->total << " in total\n\n";

    std::cout << "Current src range position: (z = " << modeling->zpos << 
                                            ", x = " << modeling->xpos << 
                                            ", y = " << modeling->ypos << ") m\n\n";
    
    std::cout << "Current rec range position: (z = " << zpos << 
                                            ", x = " << xpos << 
                                            ", y = " << ypos << ") m\n\n";

    std::cout << modeling->current_operation << "\n";
}

__global__ void image_domain_adjoint_kernel(const float * __restrict__ d_S, const float * __restrict__ d_T_src, const float * __restrict__ d_T_rec, const float * __restrict__ data,
                                            float * __restrict__ model, const float * cs_xsrc, const float * cs_ysrc, const float * cs_xrec, const float * cs_yrec, const float aperture,
                                            const float max_offset, const int sm_z, const int sm_x, const int sm_y, const float cubic_dh, const float dh, const float dt, const int cubic_nxx,
                                            const int cubic_nyy, const int cubic_nzz, const int cubic_nb, const int nx, const int ny, const int nz, const int nt, const int nsy, const int nrx)
{
    extern __shared__ float raw_smem[];
    const int total_sm_nodes = sm_z * sm_x * sm_y;

    float * sm_S  = raw_smem;
    float * sm_Ts = raw_smem + total_sm_nodes;
    float * sm_Tr = raw_smem + 2 * total_sm_nodes;

    __shared__ int cz_block_origin, cx_block_origin, cy_block_origin;

    const int tid = threadIdx.z * blockDim.y * blockDim.x + 
                    threadIdx.y * blockDim.x + 
                    threadIdx.x;
    
    const int num_threads = blockDim.x * blockDim.y * blockDim.z;

    const float inv_cubic_dh = 1.0f / cubic_dh;
    const float scale_const = 1.0f / (2.0f * M_PI);
    const int cubic_volsize = cubic_nxx * cubic_nyy * cubic_nzz;

    if (tid == 0) 
    {
        float z_phys_start = (blockIdx.x * blockDim.x) * dh;
        float x_phys_start = (blockIdx.y * blockDim.y) * dh;
        float y_phys_start = (blockIdx.z * blockDim.z) * dh;

        cz_block_origin = (int)floorf(z_phys_start * inv_cubic_dh) + cubic_nb;
        cx_block_origin = (int)floorf(x_phys_start * inv_cubic_dh) + cubic_nb;
        cy_block_origin = (int)floorf(y_phys_start * inv_cubic_dh) + cubic_nb;
    }

    __syncthreads();

    for (int idx = tid; idx < total_sm_nodes; idx += num_threads) 
    {
        int local_cz = idx % sm_z;
        int rem      = idx / sm_z;
        int local_cx = rem % sm_x;
        int local_cy = rem / sm_x;

        int global_cz = max(0, min(cz_block_origin + local_cz - 1, cubic_nzz - 1));
        int global_cx = max(0, min(cx_block_origin + local_cx - 1, cubic_nxx - 1));
        int global_cy = max(0, min(cy_block_origin + local_cy - 1, cubic_nyy - 1));

        int g_idx = global_cy * cubic_nxx * cubic_nzz + global_cx * cubic_nzz + global_cz;
        int s_idx = local_cy * sm_x * sm_z + local_cx * sm_z + local_cz;

        sm_S[s_idx] = d_S[g_idx];
    }

    __syncthreads();

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    bool valid_voxel = (i < nz) && (j < nx) && (k < ny);

    float uz = 0.0f, ux = 0.0f, uy = 0.0f;
    int sm_base_z = 0, sm_base_x = 0, sm_base_y = 0;
    float z = 0.0f, x = 0.0f, y = 0.0f;

    if (valid_voxel) 
    {
        z = __int2float_rd(i) * dh;
        x = __int2float_rd(j) * dh;
        y = __int2float_rd(k) * dh;

        int ic = (int)floorf(z * inv_cubic_dh);
        int jc = (int)floorf(x * inv_cubic_dh);
        int kc = (int)floorf(y * inv_cubic_dh);

        uz = fminf(fmaxf((z - (float)(ic) * cubic_dh) * inv_cubic_dh, 0.0f), 1.0f);
        ux = fminf(fmaxf((x - (float)(jc) * cubic_dh) * inv_cubic_dh, 0.0f), 1.0f);
        uy = fminf(fmaxf((y - (float)(kc) * cubic_dh) * inv_cubic_dh, 0.0f), 1.0f);

        sm_base_z = (ic + cubic_nb - cz_block_origin) + 1;
        sm_base_x = (jc + cubic_nb - cx_block_origin) + 1;
        sm_base_y = (kc + cubic_nb - cy_block_origin) + 1;
    }

    for (int src_csIdy = 0; src_csIdy < nsy; src_csIdy++)
    {
        int src_base = src_csIdy * cubic_volsize;

        for (int idx = tid; idx < total_sm_nodes; idx += num_threads) 
        {
            int local_cz = idx % sm_z;
            int rem      = idx / sm_z;
            int local_cx = rem % sm_x;
            int local_cy = rem / sm_x;

            int global_cz = max(0, min(cz_block_origin + local_cz - 1, cubic_nzz - 1));
            int global_cx = max(0, min(cx_block_origin + local_cx - 1, cubic_nxx - 1));
            int global_cy = max(0, min(cy_block_origin + local_cy - 1, cubic_nyy - 1));

            int g_idx = global_cy * cubic_nxx * cubic_nzz + global_cx * cubic_nzz + global_cz;
            int s_idx = local_cy * sm_x * sm_z + local_cx * sm_z + local_cz;

            sm_Ts[s_idx] = d_T_src[g_idx + src_base];
        }

        __syncthreads();

        TTD src = {0};
        if (valid_voxel) 
        {
            src = compute_TTDs(sm_Ts, sm_base_z, sm_base_x, sm_base_y, 
                              sm_z, sm_x, uz, ux, uy, cubic_dh);
        }

        for (int rec_csIdx = 0; rec_csIdx < nrx; rec_csIdx++)
        {
            int rec_base = rec_csIdx * cubic_volsize;

            for (int idx = tid; idx < total_sm_nodes; idx += num_threads) 
            {
                int local_cz = idx % sm_z;
                int rem      = idx / sm_z;
                int local_cx = rem % sm_x;
                int local_cy = rem / sm_x;

                int global_cz = max(0, min(cz_block_origin + local_cz - 1, cubic_nzz - 1));
                int global_cx = max(0, min(cx_block_origin + local_cx - 1, cubic_nxx - 1));
                int global_cy = max(0, min(cy_block_origin + local_cy - 1, cubic_nyy - 1));

                int g_idx = global_cy * cubic_nxx * cubic_nzz + global_cx * cubic_nzz + global_cz;
                int s_idx = local_cy * sm_x * sm_z + local_cx * sm_z + local_cz;

                sm_Tr[s_idx] = d_T_rec[g_idx + rec_base];
            }

            __syncthreads();

            if (valid_voxel) 
            {
                float sx = cs_xsrc[src_csIdy];
                float sy = cs_ysrc[src_csIdy];
                float rx = cs_xrec[rec_csIdx];
                float ry = cs_yrec[rec_csIdx];

                float cmpx = 0.5f * (sx + rx);
                float cmpy = 0.5f * (sy + ry);

                float offset = sqrtf((sx - rx)*(sx - rx) + (sy - ry)*(sy - ry));

                if (offset <= max_offset) 
                {
                    TTD rec = compute_TTDs(sm_Tr, sm_base_z, sm_base_x, sm_base_y, sm_z, sm_x, uz, ux, uy, cubic_dh);

                    float T = src.T + rec.T;
                    int tId = __float2int_rd(T / dt);

                    if (tId >= 0 && tId < nt) 
                    {
                        float S = tricubic(sm_S, sm_base_z, sm_base_x, sm_base_y, sm_z, sm_x, uz, ux, uy);

                        // Slowness vectors
                        float ps_x = src.dT_dx, ps_y = src.dT_dy, ps_z = src.dT_dz;
                        float pr_x = rec.dT_dx, pr_y = rec.dT_dy, pr_z = rec.dT_dz;

                        float norm_Ts = sqrtf(ps_x*ps_x + ps_y*ps_y + ps_z*ps_z) + EPS;
                        float norm_Tr = sqrtf(pr_x*pr_x + pr_y*pr_y + pr_z*pr_z) + EPS;

                        // Unit propagation direction vectors
                        float ux_s = ps_x / norm_Ts, uy_s = ps_y / norm_Ts, uz_s = ps_z / norm_Ts;
                        float ux_r = pr_x / norm_Tr, uy_r = pr_y / norm_Tr, uz_r = pr_z / norm_Tr;

                        // Obliquity factor
                        float nx_norm = 0.0f, ny_norm = 0.0f, nz_norm = -1.0f;
                        float cos_s = fminf(1.0f, fmaxf(0.001f, fabsf(ux_s*nx_norm + uy_s*ny_norm + uz_s*nz_norm)));
                        float cos_r = fminf(1.0f, fmaxf(0.001f, fabsf(ux_r*nx_norm + uy_r*ny_norm + uz_r*nz_norm)));
                        float obliquity = 0.5f * sqrtf(cos_s + cos_r);

                        // Stable Beylkin Opening Angle Jacobian
                        float dot_sr = ux_s*ux_r + uy_s*uy_r + uz_s*uz_r;
                        float cos_psi = fminf(1.0f, fmaxf(-1.0f, dot_sr));
                        float sin_psi = sqrtf(1.0f - cos_psi * cos_psi);

                        float sum_p_mag = sqrtf((ps_x + pr_x)*(ps_x + pr_x) + 
                                                (ps_y + pr_y)*(ps_y + pr_y) + 
                                                (ps_z + pr_z)*(ps_z + pr_z));

                        float Jterm = sum_p_mag * sin_psi;

                        // Geometrical Spreading
                        float R_s = fmaxf(src.T / S, EPS);
                        float R_r = fmaxf(rec.T / S, EPS);
                        float G = 1.0f / sqrtf(R_s * R_r);

                        // Reflection Angle Weight
                        float theta_s = acosf(fminf(1.0f, fmaxf(-1.0f, ux_s*nx_norm + uy_s*ny_norm + uz_s*nz_norm)));
                        float theta_r = acosf(fminf(1.0f, fmaxf(-1.0f, ux_r*nx_norm + uy_r*ny_norm + uz_r*nz_norm)));
                        float theta = 0.5f * (theta_s + theta_r);
                        float R = 1.0f + 0.2f * cosf(theta);
                        
                        // Aperture Gaussian Tapering
                        float sigma = tanf(aperture * M_PI / 180.0f) * z;
                        float par_x = ((x - cmpx) / (sigma + EPS)) * ((x - cmpx) / (sigma + EPS));
                        float par_y = ((y - cmpy) / (sigma + EPS)) * ((y - cmpy) / (sigma + EPS));
                        float taper = expf(-0.5f * (par_x + par_y));

                        // Total Kernel Weights
                        float weights = scale_const * taper * G * R * obliquity * Jterm;

                        int mId = i + j * nz + k * nx * nz;

                        int traceId = src_csIdy * nrx + rec_csIdx;
                        int tId_base = traceId * nt;
                
                        atomicAdd(&model[mId], weights * data[tId + tId_base]);
                    }
                }
            }

            __syncthreads();
        }

        __syncthreads();
    }
}

__global__ void image_domain_forward_kernel(const float * __restrict__ d_S, const float * __restrict__ d_T_src, const float * __restrict__ d_T_rec, float * __restrict__ data,
                                            const float * __restrict__ model, const float * cs_xsrc, const float * cs_ysrc, const float * cs_xrec, const float * cs_yrec, const float aperture,
                                            const float max_offset, const int sm_z, const int sm_x, const int sm_y, const float cubic_dh, const float dh, const float dt, const int cubic_nxx,
                                            const int cubic_nyy, const int cubic_nzz, const int cubic_nb, const int nx, const int ny, const int nz, const int nt, const int nsy, const int nrx)
{
    extern __shared__ float raw_smem[];
    const int total_sm_nodes = sm_z * sm_x * sm_y;

    float * sm_S  = raw_smem;
    float * sm_Ts = raw_smem + total_sm_nodes;
    float * sm_Tr = raw_smem + 2 * total_sm_nodes;

    __shared__ int cz_block_origin, cx_block_origin, cy_block_origin;

    const int tid = threadIdx.z * blockDim.y * blockDim.x + 
                    threadIdx.y * blockDim.x + 
                    threadIdx.x;
    
    const int num_threads = blockDim.x * blockDim.y * blockDim.z;

    const float inv_cubic_dh = 1.0f / cubic_dh;
    const float scale_const = 1.0f / (2.0f * M_PI);
    const int cubic_volsize = cubic_nxx * cubic_nyy * cubic_nzz;

    if (tid == 0) 
    {
        float z_phys_start = (blockIdx.x * blockDim.x) * dh;
        float x_phys_start = (blockIdx.y * blockDim.y) * dh;
        float y_phys_start = (blockIdx.z * blockDim.z) * dh;

        cz_block_origin = (int)floorf(z_phys_start * inv_cubic_dh) + cubic_nb;
        cx_block_origin = (int)floorf(x_phys_start * inv_cubic_dh) + cubic_nb;
        cy_block_origin = (int)floorf(y_phys_start * inv_cubic_dh) + cubic_nb;
    }

    __syncthreads();

    for (int idx = tid; idx < total_sm_nodes; idx += num_threads) 
    {
        int local_cz = idx % sm_z;
        int rem      = idx / sm_z;
        int local_cx = rem % sm_x;
        int local_cy = rem / sm_x;

        int global_cz = max(0, min(cz_block_origin + local_cz - 1, cubic_nzz - 1));
        int global_cx = max(0, min(cx_block_origin + local_cx - 1, cubic_nxx - 1));
        int global_cy = max(0, min(cy_block_origin + local_cy - 1, cubic_nyy - 1));

        int g_idx = global_cy * cubic_nxx * cubic_nzz + global_cx * cubic_nzz + global_cz;
        int s_idx = local_cy * sm_x * sm_z + local_cx * sm_z + local_cz;

        sm_S[s_idx] = d_S[g_idx];
    }

    __syncthreads();

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    int k = blockIdx.z * blockDim.z + threadIdx.z;

    bool valid_voxel = (i < nz) && (j < nx) && (k < ny);

    float uz = 0.0f, ux = 0.0f, uy = 0.0f;
    int sm_base_z = 0, sm_base_x = 0, sm_base_y = 0;
    float z = 0.0f, x = 0.0f, y = 0.0f;

    if (valid_voxel) 
    {
        z = __int2float_rd(i) * dh;
        x = __int2float_rd(j) * dh;
        y = __int2float_rd(k) * dh;

        int ic = (int)floorf(z * inv_cubic_dh);
        int jc = (int)floorf(x * inv_cubic_dh);
        int kc = (int)floorf(y * inv_cubic_dh);

        uz = fminf(fmaxf((z - (float)(ic) * cubic_dh) * inv_cubic_dh, 0.0f), 1.0f);
        ux = fminf(fmaxf((x - (float)(jc) * cubic_dh) * inv_cubic_dh, 0.0f), 1.0f);
        uy = fminf(fmaxf((y - (float)(kc) * cubic_dh) * inv_cubic_dh, 0.0f), 1.0f);

        sm_base_z = (ic + cubic_nb - cz_block_origin) + 1;
        sm_base_x = (jc + cubic_nb - cx_block_origin) + 1;
        sm_base_y = (kc + cubic_nb - cy_block_origin) + 1;
    }

    for (int src_csIdy = 0; src_csIdy < nsy; src_csIdy++)
    {
        int src_base = src_csIdy * cubic_volsize;

        for (int idx = tid; idx < total_sm_nodes; idx += num_threads) 
        {
            int local_cz = idx % sm_z;
            int rem      = idx / sm_z;
            int local_cx = rem % sm_x;
            int local_cy = rem / sm_x;

            int global_cz = max(0, min(cz_block_origin + local_cz - 1, cubic_nzz - 1));
            int global_cx = max(0, min(cx_block_origin + local_cx - 1, cubic_nxx - 1));
            int global_cy = max(0, min(cy_block_origin + local_cy - 1, cubic_nyy - 1));

            int g_idx = global_cy * cubic_nxx * cubic_nzz + global_cx * cubic_nzz + global_cz;
            int s_idx = local_cy * sm_x * sm_z + local_cx * sm_z + local_cz;

            sm_Ts[s_idx] = d_T_src[g_idx + src_base];
        }

        __syncthreads();

        TTD src = {0};
        if (valid_voxel) 
            src = compute_TTDs(sm_Ts, sm_base_z, sm_base_x, sm_base_y, sm_z, sm_x, uz, ux, uy, cubic_dh);

        for (int rec_csIdx = 0; rec_csIdx < nrx; rec_csIdx++)
        {
            int rec_base = rec_csIdx * cubic_volsize;

            for (int idx = tid; idx < total_sm_nodes; idx += num_threads) 
            {
                int local_cz = idx % sm_z;
                int rem      = idx / sm_z;
                int local_cx = rem % sm_x;
                int local_cy = rem / sm_x;

                int global_cz = max(0, min(cz_block_origin + local_cz - 1, cubic_nzz - 1));
                int global_cx = max(0, min(cx_block_origin + local_cx - 1, cubic_nxx - 1));
                int global_cy = max(0, min(cy_block_origin + local_cy - 1, cubic_nyy - 1));

                int g_idx = global_cy * cubic_nxx * cubic_nzz + global_cx * cubic_nzz + global_cz;
                int s_idx = local_cy * sm_x * sm_z + local_cx * sm_z + local_cz;

                sm_Tr[s_idx] = d_T_rec[g_idx + rec_base];
            }

            __syncthreads();

            if (valid_voxel) 
            {
                float sx = cs_xsrc[src_csIdy];
                float sy = cs_ysrc[src_csIdy];
                float rx = cs_xrec[rec_csIdx];
                float ry = cs_yrec[rec_csIdx];

                float cmpx = 0.5f * (sx + rx);
                float cmpy = 0.5f * (sy + ry);

                float offset = sqrtf((sx - rx)*(sx - rx) + (sy - ry)*(sy - ry));

                if (offset <= max_offset) 
                {
                    TTD rec = compute_TTDs(sm_Tr, sm_base_z, sm_base_x, sm_base_y, sm_z, sm_x, uz, ux, uy, cubic_dh);

                    float T = src.T + rec.T;
                    int tId = __float2int_rd(T / dt);

                    if (tId >= 0 && tId < nt) 
                    {
                        float S = tricubic(sm_S, sm_base_z, sm_base_x, sm_base_y, sm_z, sm_x, uz, ux, uy);

                        // Slowness vectors
                        float ps_x = src.dT_dx, ps_y = src.dT_dy, ps_z = src.dT_dz;
                        float pr_x = rec.dT_dx, pr_y = rec.dT_dy, pr_z = rec.dT_dz;

                        float norm_Ts = sqrtf(ps_x*ps_x + ps_y*ps_y + ps_z*ps_z) + EPS;
                        float norm_Tr = sqrtf(pr_x*pr_x + pr_y*pr_y + pr_z*pr_z) + EPS;

                        // Unit propagation direction vectors
                        float ux_s = ps_x / norm_Ts, uy_s = ps_y / norm_Ts, uz_s = ps_z / norm_Ts;
                        float ux_r = pr_x / norm_Tr, uy_r = pr_y / norm_Tr, uz_r = pr_z / norm_Tr;

                        // Obliquity factor
                        float nx_norm = 0.0f, ny_norm = 0.0f, nz_norm = -1.0f;
                        float cos_s = fminf(1.0f, fmaxf(0.001f, fabsf(ux_s*nx_norm + uy_s*ny_norm + uz_s*nz_norm)));
                        float cos_r = fminf(1.0f, fmaxf(0.001f, fabsf(ux_r*nx_norm + uy_r*ny_norm + uz_r*nz_norm)));
                        float obliquity = 0.5f * sqrtf(cos_s + cos_r);

                        // Stable Beylkin Opening Angle Jacobian
                        float dot_sr = ux_s*ux_r + uy_s*uy_r + uz_s*uz_r;
                        float cos_psi = fminf(1.0f, fmaxf(-1.0f, dot_sr));
                        float sin_psi = sqrtf(1.0f - cos_psi * cos_psi);

                        float sum_p_mag = sqrtf((ps_x + pr_x)*(ps_x + pr_x) + 
                                                (ps_y + pr_y)*(ps_y + pr_y) + 
                                                (ps_z + pr_z)*(ps_z + pr_z));

                        float Jterm = sum_p_mag * sin_psi;

                        // Geometrical Spreading
                        float R_s = fmaxf(src.T / S, EPS);
                        float R_r = fmaxf(rec.T / S, EPS);
                        float G = 1.0f / sqrtf(R_s * R_r);

                        // Reflection Angle Weight
                        float theta_s = acosf(fminf(1.0f, fmaxf(-1.0f, ux_s*nx_norm + uy_s*ny_norm + uz_s*nz_norm)));
                        float theta_r = acosf(fminf(1.0f, fmaxf(-1.0f, ux_r*nx_norm + uy_r*ny_norm + uz_r*nz_norm)));
                        float theta = 0.5f * (theta_s + theta_r);
                        float R = 1.0f + 0.2f * cosf(theta);
                        
                        // Aperture Gaussian Tapering
                        float sigma = tanf(aperture * M_PI / 180.0f) * z;
                        float par_x = ((x - cmpx) / (sigma + EPS)) * ((x - cmpx) / (sigma + EPS));
                        float par_y = ((y - cmpy) / (sigma + EPS)) * ((y - cmpy) / (sigma + EPS));
                        float taper = expf(-0.5f * (par_x + par_y));

                        // Total Kernel Weights
                        float weights = scale_const * taper * G * R * obliquity * Jterm;

                        int mId = i + j * nz + k * nx * nz;

                        int traceId = src_csIdy * nrx + rec_csIdx;
                        int tId_base = traceId * nt;
                
                        atomicAdd(&data[tId + tId_base], weights * model[mId]);
                    }
                }
            }

            __syncthreads();
        }

        __syncthreads();
    }
}

__global__ void angle_domain_adjoint_kernel()
{


}

__global__ void angle_domain_forward_kernel()
{


}

__device__ __forceinline__ void get_catmull_weights(float u, float w[4], float dw[4], float d2w[4]) 
{
    float u2 = u * u;
    float u3 = u2 * u;

    w[0] = -0.5f * u3 + u2 - 0.5f * u;
    w[1] =  1.5f * u3 - 2.5f * u2 + 1.0f;
    w[2] = -1.5f * u3 + 2.0f * u2 + 0.5f * u;
    w[3] =  0.5f * u3 - 0.5f * u2;

    dw[0] = -1.5f * u2 + 2.0f * u - 0.5f;
    dw[1] =  4.5f * u2 - 5.0f * u;
    dw[2] = -4.5f * u2 + 4.0f * u + 0.5f;
    dw[3] =  1.5f * u2 - u;

    d2w[0] = -3.0f * u + 2.0f;
    d2w[1] =  9.0f * u - 5.0f;
    d2w[2] = -9.0f * u + 4.0f;
    d2w[3] =  3.0f * u - 1.0f;
}

__device__ __forceinline__ TTD compute_TTDs(const float * __restrict__ smem_buffer, int base_z, int base_x, int base_y, 
                                            int smem_z, int smem_x, float uz, float ux, float uy, float cubic_dh) 
{
    float wz[4], dwz[4], d2wz[4];
    float wx[4], dwx[4], d2wx[4];
    float wy[4], dwy[4], d2wy[4];

    get_catmull_weights(uz, wz, dwz, d2wz);
    get_catmull_weights(ux, wx, dwx, d2wx);
    get_catmull_weights(uy, wy, dwy, d2wy);

    const int stride_x = smem_z;
    const int stride_y = smem_x * smem_z;

    float inv_z  = 1.0f / cubic_dh;
    float inv_x  = 1.0f / cubic_dh;
    float inv_y  = 1.0f / cubic_dh;

    float inv_z2 = inv_z * inv_z;
    float inv_x2 = inv_x * inv_x;
    float inv_y2 = inv_y * inv_y;

    float inv_xz = inv_x * inv_z;
    float inv_yz = inv_y * inv_z;
    float inv_xy = inv_x * inv_y;

    float R0[4][4], R1[4][4], R2[4][4];

    #pragma unroll
    for (int dy = 0; dy < 4; ++dy) 
    {
        int y_offset = (base_y + dy - 1) * stride_y;

        #pragma unroll
        for (int dx = 0; dx < 4; ++dx) 
        {
            int offset = y_offset + (base_x + dx - 1) * stride_x + (base_z - 1);

            float r0 = 0.0f, r1 = 0.0f, r2 = 0.0f;
        
            #pragma unroll
            for (int dz = 0; dz < 4; ++dz) 
            {
                float val = smem_buffer[offset + dz];
                
                r0 += val * wz[dz];   
                r1 += val * dwz[dz];  
                r2 += val * d2wz[dz]; 
            }

            R0[dy][dx] = r0;
            R1[dy][dx] = r1;
            R2[dy][dx] = r2;
        }
    }

    float RX00[4], RX01[4], RX02[4];
    float RX10[4], RX11[4];
    float RX20[4];

    #pragma unroll
    for (int dy = 0; dy < 4; ++dy) 
    {
        float rx00 = 0.0f, rx01 = 0.0f, rx02 = 0.0f;
        float rx10 = 0.0f, rx11 = 0.0f;
        float rx20 = 0.0f;

        #pragma unroll
        for (int dx = 0; dx < 4; ++dx) 
        {
            float r0 = R0[dy][dx];
            float r1 = R1[dy][dx];
            float r2 = R2[dy][dx];

            float w = wx[dx], dw = dwx[dx], d2w = d2wx[dx];

            rx00 += r0 * w;   
            rx01 += r0 * dw;  
            rx02 += r0 * d2w; 

            rx10 += r1 * w;   
            rx11 += r1 * dw;  

            rx20 += r2 * w;   
        }

        RX00[dy] = rx00; RX01[dy] = rx01; RX02[dy] = rx02;
        RX10[dy] = rx10; RX11[dy] = rx11;
        RX20[dy] = rx20;
    }

    TTD ttd = {0};

    #pragma unroll
    for (int dy = 0; dy < 4; ++dy) 
    {
        float w = wy[dy], dw = dwy[dy], d2w = d2wy[dy];

        ttd.T        += RX00[dy] * w;

        ttd.dT_dz    += RX10[dy] * w;
        ttd.dT_dx    += RX01[dy] * w;
        ttd.dT_dy    += RX00[dy] * dw;

        ttd.d2T_dz2  += RX20[dy] * w;
        ttd.d2T_dx2  += RX02[dy] * w;
        ttd.d2T_dy2  += RX00[dy] * d2w;

        ttd.d2T_dxdz += RX11[dy] * w;
        ttd.d2T_dydz += RX10[dy] * dw;
        ttd.d2T_dxdy += RX01[dy] * dw;
    }

    ttd.dT_dz    *= inv_z;
    ttd.dT_dx    *= inv_x;
    ttd.dT_dy    *= inv_y;

    ttd.d2T_dz2  *= inv_z2;
    ttd.d2T_dx2  *= inv_x2;
    ttd.d2T_dy2  *= inv_y2;

    ttd.d2T_dxdz *= inv_xz;
    ttd.d2T_dydz *= inv_yz;
    ttd.d2T_dxdy *= inv_xy;

    return ttd;
}

__device__ __forceinline__ float cubic1d(float p0, float p1, float p2, float p3, float u) 
{
    float c0 = p1;
    float c1 = 0.5f * (p2 - p0);
    float c2 = p0 - 2.5f * p1 + 2.0f * p2 - 0.5f * p3;
    float c3 = 0.5f * (p3 - p0) + 1.5f * (p1 - p2);
    
    return ((c3 * u + c2) * u + c1) * u + c0;
}

__device__ __forceinline__ float tricubic(const float * __restrict__ smem_buffer, int base_z, int base_x, int base_y, 
                                          int smem_z, int smem_x, float uz, float ux, float uy) 
{
    const int stride_x = smem_z;
    const int stride_y = smem_x * smem_z;

    float val_y[4];

    #pragma unroll
    for (int dy = 0; dy < 4; ++dy) 
    {
        int y_offset = (base_y + dy - 1) * stride_y;
        float val_x[4]; 

        #pragma unroll
        for (int dx = 0; dx < 4; ++dx) 
        {
            int xy_offset = y_offset + (base_x + dx - 1) * stride_x + (base_z - 1);

            float p0 = smem_buffer[xy_offset];
            float p1 = smem_buffer[xy_offset + 1];
            float p2 = smem_buffer[xy_offset + 2];
            float p3 = smem_buffer[xy_offset + 3];

            val_x[dx] = cubic1d(p0, p1, p2, p3, uz);
        }

        val_y[dy] = cubic1d(val_x[0], val_x[1], val_x[2], val_x[3], ux);
    }

    return cubic1d(val_y[0], val_y[1], val_y[2], val_y[3], uy);
}


