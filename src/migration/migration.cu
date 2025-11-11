# include "migration.cuh"

void Migration::set_parameters()
{
    nt = std::stoi(catch_parameter("time_samples", parameters));
    dt = std::stof(catch_parameter("time_spacing", parameters));
    da = std::stof(catch_parameter("mig_angle_spacing", parameters));

    max_angle = std::stof(catch_parameter("mig_max_angle", parameters));
    max_it = std::stoi(catch_parameter("mig_max_iteration", parameters));
    residuo_folder = catch_parameter("mig_residuo_folder", parameters);

    fmax = std::stof(catch_parameter("max_frequency", parameters));

    input_data_folder = catch_parameter("mig_data_folder", parameters);
    input_data_prefix = catch_parameter("mig_data_prefix", parameters);
    
    tables_folder = catch_parameter("mig_tables_folder", parameters);
    images_folder = catch_parameter("mig_images_folder", parameters);
    gathers_folder = catch_parameter("mig_gathers_folder", parameters);

    anisotropy = str2bool(catch_parameter("anisotropy", parameters));

    if (anisotropy) modeling = new Eikonal_ANI(); 
               else modeling = new Eikonal_ISO();

    modeling->parameters = parameters;
    modeling->set_parameters();

    nang = (int)(max_angle / da) + 1;

    nBlocks = (int)((modeling->volsize + NTHREADS - 1) / NTHREADS);

    set_wavelet();
    set_gathers();
    
    set_migration();

    d_samples = nt*modeling->geometry->nrec*modeling->geometry->nsrc;

    h_data = new float[nt]();    
    h_model = new float[m_samples]();

    h_Ts = new float[modeling->volsize]();
    h_Tr = new float[modeling->volsize]();
    
    seismic = new float[nt*modeling->geometry->nrec]();

    cudaMalloc((void**)&(d_data), nt*sizeof(float));    
    cudaMalloc((void**)&(d_model), m_samples*sizeof(float));
    
    cudaMalloc((void**)&(d_Ts), modeling->volsize*sizeof(float));
    cudaMalloc((void**)&(d_Tr), modeling->volsize*sizeof(float));
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

void Migration::set_gathers()
{
    ds = modeling->geometry->ysrc[1] - modeling->geometry->ysrc[0];
    dr = modeling->geometry->yrec[1] - modeling->geometry->yrec[0];

    float dCMP = 0.5f*min(ds,dr);

    minCMPx = 1e6f;
    minCMPy = 1e6f;

    float maxCMPy =-1e6f;  
    float maxCMPx =-1e6f;  

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

void Migration::set_src_domain()
{
    keyword = "source";
    total = std::to_string(modeling->geometry->nsrc); 
}

void Migration::set_current_src()
{
    modeling->sx = modeling->geometry->xsrc[modeling->srcId];
    modeling->sy = modeling->geometry->ysrc[modeling->srcId];
    modeling->sz = modeling->geometry->zsrc[modeling->srcId];

    current = std::to_string(modeling->srcId+1);
    
    xpos = format1Decimal(modeling->sx);
    ypos = format1Decimal(modeling->sy);
    zpos = format1Decimal(modeling->sz);
}

void Migration::set_rec_domain()
{
    keyword = "receiver";
    total = std::to_string(modeling->geometry->nrec); 
}

void Migration::set_current_rec()
{
    modeling->sx = modeling->geometry->xrec[modeling->recId];
    modeling->sy = modeling->geometry->yrec[modeling->recId];
    modeling->sz = modeling->geometry->zrec[modeling->recId];

    current = std::to_string(modeling->recId+1);
    
    xpos = format1Decimal(modeling->sx);
    ypos = format1Decimal(modeling->sy);
    zpos = format1Decimal(modeling->sz);
}

void Migration::set_src_travel_times()
{
    set_src_domain();

    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nsrc; modeling->srcId++)
    {
        set_current_src();
    
        current_operation = "Computing " + keyword + " travel time matrices";

        show_information();

        modeling->time_propagation();
        
        cudaMemcpy(h_Ts, modeling->d_T, modeling->volsize*sizeof(float), cudaMemcpyDeviceToHost);

        export_binary_float(tables_folder + "eikonal_src_" + std::to_string(modeling->srcId+1) + ".bin", h_Ts, modeling->volsize);
    }
}

void Migration::set_rec_travel_times()
{
    set_rec_domain();
    
    for (modeling->recId = 0; modeling->recId < modeling->geometry->nrec; modeling->recId++)
    {
        set_current_rec();

        current_operation = "Computing " + keyword + " travel time matrices";

        show_information();

        modeling->time_propagation();
        
        cudaMemcpy(h_Tr, modeling->d_T, modeling->volsize*sizeof(float), cudaMemcpyDeviceToHost);

        export_binary_float(tables_folder + "eikonal_rec_" + std::to_string(modeling->recId+1) + ".bin", h_Tr, modeling->volsize);
    }
}

void Migration::prepare_convolution()
{
    nfft = nextpow2(nt);

    time_trace = (double *) fftw_malloc(nfft*sizeof(double));
    time_wavelet = (double *) fftw_malloc(nfft*sizeof(double));

    freq_trace = (fftw_complex *) fftw_malloc(nfft*sizeof(fftw_complex));
    freq_wavelet = (fftw_complex *) fftw_malloc(nfft*sizeof(fftw_complex));

    trace_forward_plan = fftw_plan_dft_r2c_1d(nt, time_trace, freq_trace, FFTW_ESTIMATE);
    trace_inverse_plan = fftw_plan_dft_c2r_1d(nt, freq_trace, time_trace, FFTW_ESTIMATE);
    wavelet_forward_plan = fftw_plan_dft_r2c_1d(nt, time_wavelet, freq_wavelet, FFTW_ESTIMATE);

    for (int tId = 0; tId < nfft; tId++)
    {
        time_trace[tId] = 0.0;
        time_wavelet[tId] = (tId < nw) ? wavelet[tId] : 0.0;
    }

    fftw_execute(wavelet_forward_plan);
}

void Migration::show_information()
{
    auto clear = system("clear");
    
    std::cout << "-------------------------------------------------------------------\n";
    std::cout << " \033[34mSeisFAT3D\033[0;0m --------------------------------------------------------\n";
    std::cout << "-------------------------------------------------------------------\n\n";

    std::cout << "Model dimensions: (z = " << (modeling->nz - 1)*modeling->dz << 
                                  ", x = " << (modeling->nx - 1)*modeling->dx << 
                                  ", y = " << (modeling->ny - 1)*modeling->dy << ") m\n\n";

    std::cout << "Running " << keyword << " " << current << " of " << total << " in total\n\n";

    std::cout << "Current " << keyword << " position: (z = " << zpos << 
                                                    ", x = " << xpos << 
                                                    ", y = " << ypos << ") m\n\n";
    
    std::cout << current_operation << "\n";                                                   
}

void Migration::adjoint_convolution()
{
    for (int tId = 0; tId < nt; tId++)
    {
        time_trace[tId] = (double)h_data[tId];
        h_data[tId] = 0.0f;
    }
    
    fftw_execute(trace_forward_plan);

    for (int fId = 0; fId < nfft; fId++)
    {
        double a_re = freq_trace[fId][0];
        double a_im = freq_trace[fId][1];

        double b_re = freq_wavelet[fId][0];
        double b_im = freq_wavelet[fId][1];

        freq_trace[fId][0] = a_re * b_re + a_im * b_im;  
        freq_trace[fId][1] = a_im * b_re - a_re * b_im;  
    }

    fftw_execute(trace_inverse_plan);

    for (int tId = nw/2 + nw/10; tId < nt; tId++)
        h_data[tId] = (float)time_trace[tId - nw/2 - nw/10] / nfft;
}

void Migration::forward_convolution()
{
    for (int tId = 0; tId < nt; tId++)
    {
        time_trace[tId] = (double)h_data[tId];
        h_data[tId] = 0.0f;
    }

    fftw_execute(trace_forward_plan);

    for (int fId = 0; fId < nfft; fId++)
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
        h_data[tId] = (float)time_trace[tId + nw/2 + nw/10] / nfft / nfft;
}

void Migration::dot_product_test()
{
    set_src_travel_times();
    set_rec_travel_times();
    prepare_convolution();
    
    auto trash = system("clear");
    std::cout << "Computing dot product test!\n\n";

    // d1 = new float[d_samples]();
    // d2 = new float[d_samples]();

    // m1 = new float[m_samples]();
    // m2 = new float[m_samples]();

    // int minValue =-100;
    // int maxValue = 100;

    // std::mt19937 prng(std::random_device{}());
    // std::uniform_int_distribution<int> dist(minValue, maxValue);

    // for (int mId = 0; mId < m_samples; mId++)
    //     m1[mId] = dist(prng);

    // for (int dId = 0; dId < d_samples; dId++)
    //     d2[dId] = dist(prng);
    
    // cudaMemcpy(d_model, m1, m_samples*sizeof(float), cudaMemcpyHostToDevice);

    // for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nsrc; modeling->srcId++)
    // {     
    //     import_binary_float(tables_folder + "eikonal_src_" + std::to_string(modeling->srcId+1) + ".bin", h_Ts, modeling->matsize);
    //     cudaMemcpy(d_Ts, h_Ts, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);
        
    //     int spreadId = 0;
        
    //     for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
    //     {
    //         cmpId = spreadId + 2.0f*(ds/dr)*modeling->srcId;                              
    //         cmp = 0.5f*(modeling->geometry->xsrc[modeling->srcId] + 
    //                     modeling->geometry->xrec[modeling->recId]);

    //         import_binary_float(tables_folder + "eikonal_rec_" + std::to_string(modeling->recId+1) + ".bin", h_Tr, modeling->matsize);
    //         cudaMemcpy(d_Tr, h_Tr, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

    //         cudaMemset(d_data, 0.0f, nt * sizeof(float));

    //         perform_forward();

    //         cudaMemcpy(h_data, d_data, nt * sizeof(float), cudaMemcpyDeviceToHost);

    //         forward_convolution();

    //         for (int tId = 0; tId < nt; tId++)
    //         {
    //             int index = tId + spreadId*nt + modeling->srcId*modeling->max_spread*nt;
                
    //             d1[index] = h_data[tId];    
    //         }                
        
    //         ++spreadId;
    //     }
    // }

    // cudaMemset(d_model, 0.0f, m_samples*sizeof(float));

    // for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nsrc; modeling->srcId++)
    // { 
    //     modeling->sx = modeling->geometry->xsrc[modeling->srcId];

    //     import_binary_float(tables_folder + "eikonal_src_" + std::to_string(modeling->srcId+1) + ".bin", h_Ts, modeling->matsize);
    //     cudaMemcpy(d_Ts, h_Ts, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

    //     int spreadId = 0;
        
    //     for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
    //     {
    //         cmp = 0.5f*(modeling->sx + modeling->geometry->xrec[modeling->recId]);
    //         cmpId = spreadId + 2.0f*(ds/dr)*modeling->srcId;                              

    //         import_binary_float(tables_folder + "eikonal_rec_" + std::to_string(modeling->recId+1) + ".bin", h_Tr, modeling->matsize);
    //         cudaMemcpy(d_Tr, h_Tr, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

    //         for (int tId = 0; tId < nt; tId++)
    //             h_data[tId] = d2[tId + spreadId*nt + modeling->srcId*modeling->max_spread*nt];

    //         adjoint_convolution();

    //         cudaMemcpy(d_data, h_data, nt * sizeof(float), cudaMemcpyHostToDevice);
            
    //         perform_adjoint();

    //         ++spreadId;
    //     }
    // }

    // cudaMemcpy(m2, d_model, m_samples*sizeof(float), cudaMemcpyDeviceToHost);        

    // double dot_m = 0.0;
    // double dot_d = 0.0;

    // for (int mId = 0; mId < m_samples; mId++)
    //     dot_m += m1[mId]*m2[mId];
    
    // for (int dId = 0; dId < d_samples; dId++)
    //     dot_d += d1[dId]*d2[dId];

    // double r = (dot_d - dot_m) / (dot_d + dot_m);
    
    // std::cout << "<m1, m2> = " << dot_m << std::endl;
    // std::cout << "<d1, d2> = " << dot_d << std::endl; 
    // std::cout << "residuo = " << r << std::endl;
}

__global__ void image_domain_adjoint_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dx, float dz, float dt, int nxx, int nzz, int nt, int nb)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int i = (int)(index % nzz);
    int j = (int)(index / nzz);

    if ((i > nb) && (i < nzz-nb) && (j > nb) && (j < nxx-nb))
    {
        float T = Ts[index] + Tr[index]; 
        int tId = __float2int_rd(T / dt);

        int nz = (nzz - 2*nb);

        if (tId < nt) 
        {
            float dTs_dx = 0.5f*(Ts[i + (j+1)*nzz] - Ts[i + (j-1)*nzz]) / dx;
            float dTs_dz = 0.5f*(Ts[(i+1) + j*nzz] - Ts[(i-1) + j*nzz]) / dz;

            float d2Ts_dx2 = (Ts[i + (j+1)*nzz] - 2.0f*Ts[index] + Ts[i + (j-1)*nzz]) / (dx*dx); 
            float d2Ts_dz2 = (Ts[(i+1) + j*nzz] - 2.0f*Ts[index] + Ts[(i-1) + j*nzz]) / (dz*dz);
           
            float d2Ts_dxdz = (Ts[(i+1) + (j+1)*nzz] - Ts[(i-1) + (j+1)*nzz] - Ts[(i+1) + (j-1)*nzz] + Ts[(i-1) + (j-1)*nzz]) / (4.0f*dx*dz);

            float dTr_dx = 0.5f*(Tr[i + (j+1)*nzz] - Tr[i + (j-1)*nzz]) / dx;
            float dTr_dz = 0.5f*(Tr[(i+1) + j*nzz] - Tr[(i-1) + j*nzz]) / dz;

            float d2Tr_dx2 = (Tr[i + (j+1)*nzz] - 2.0f*Tr[index] + Tr[i + (j-1)*nzz]) / (dx*dx); 
            float d2Tr_dz2 = (Tr[(i+1) + j*nzz] - 2.0f*Tr[index] + Tr[(i-1) + j*nzz]) / (dz*dz);
           
            float d2Tr_dxdz = (Tr[(i+1) + (j+1)*nzz] - Tr[(i-1) + (j+1)*nzz] - Tr[(i+1) + (j-1)*nzz] + Tr[(i-1) + (j-1)*nzz]) / (4.0f*dx*dz);

            float norm_Ts = sqrtf(dTs_dx*dTs_dx + dTs_dz*dTs_dz) + EPS;
            float norm_Tr = sqrtf(dTr_dx*dTr_dx + dTr_dz*dTr_dz) + EPS;
            
            float ux_s = dTs_dx / norm_Ts; float uz_s = dTs_dz / norm_Ts;            
            float ux_r = dTr_dx / norm_Tr; float uz_r = dTr_dz / norm_Tr;            

            float nx_norm = 0.0f, nz_norm = -1.0f; 
            float cos_s = fabs(ux_s*nx_norm + uz_s*nz_norm);
            float cos_r = fabs(ux_r*nx_norm + uz_r*nz_norm);

            float a = d2Ts_dx2  + d2Tr_dx2;
            float b = d2Ts_dxdz + d2Tr_dxdz;
            float c = d2Ts_dz2  + d2Tr_dz2;

            float detH = a*c - b*b;

            float J = 1.0f / sqrtf(fabsf(detH) + EPS); 

            float R_s = max(Ts[index] / S[index], EPS);
            float R_r = max(Tr[index] / S[index], EPS);

            float G = 1.0f / sqrt(R_s * R_r);

            float theta = acosf(min(1.0f, max(-1.0f, ux_s*nx_norm + uz_s*nz_norm)));

            float R = 1.0f + 0.2f*cos(theta);

            float weights = 1.0f / (2.0f * M_PI) * sqrt(max(0.0f, cos_s) * max(0.0f, cos_r)) * G * J * R;            

            int mId = (i - nb) + (j - nb)*nz;
    
            atomicAdd(&model[mId], weights * data[tId]);
        }
    }    
}

__global__ void image_domain_forward_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dx, float dz, float dt, int nxx, int nzz, int nt, int nb)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int i = (int)(index % nzz);
    int j = (int)(index / nzz);

    if ((i > nb) && (i < nzz-nb) && (j > nb) && (j < nxx-nb))
    {
        float T = Ts[index] + Tr[index]; 
        int tId = __float2int_rd(T / dt);

        int nz = (nzz - 2*nb);

        if (tId < nt) 
        {
            float dTs_dx = 0.5f*(Ts[i + (j+1)*nzz] - Ts[i + (j-1)*nzz]) / dx;
            float dTs_dz = 0.5f*(Ts[(i+1) + j*nzz] - Ts[(i-1) + j*nzz]) / dz;

            float d2Ts_dx2 = (Ts[i + (j+1)*nzz] - 2.0f*Ts[index] + Ts[i + (j-1)*nzz]) / (dx*dx); 
            float d2Ts_dz2 = (Ts[(i+1) + j*nzz] - 2.0f*Ts[index] + Ts[(i-1) + j*nzz]) / (dz*dz);
           
            float d2Ts_dxdz = (Ts[(i+1) + (j+1)*nzz] - Ts[(i-1) + (j+1)*nzz] - Ts[(i+1) + (j-1)*nzz] + Ts[(i-1) + (j-1)*nzz]) / (4.0f*dx*dz);

            float dTr_dx = 0.5f*(Tr[i + (j+1)*nzz] - Tr[i + (j-1)*nzz]) / dx;
            float dTr_dz = 0.5f*(Tr[(i+1) + j*nzz] - Tr[(i-1) + j*nzz]) / dz;

            float d2Tr_dx2 = (Tr[i + (j+1)*nzz] - 2.0f*Tr[index] + Tr[i + (j-1)*nzz]) / (dx*dx); 
            float d2Tr_dz2 = (Tr[(i+1) + j*nzz] - 2.0f*Tr[index] + Tr[(i-1) + j*nzz]) / (dz*dz);
           
            float d2Tr_dxdz = (Tr[(i+1) + (j+1)*nzz] - Tr[(i-1) + (j+1)*nzz] - Tr[(i+1) + (j-1)*nzz] + Tr[(i-1) + (j-1)*nzz]) / (4.0f*dx*dz);

            float norm_Ts = sqrtf(dTs_dx*dTs_dx + dTs_dz*dTs_dz) + EPS;
            float norm_Tr = sqrtf(dTr_dx*dTr_dx + dTr_dz*dTr_dz) + EPS;
            
            float ux_s = dTs_dx / norm_Ts; float uz_s = dTs_dz / norm_Ts;            
            float ux_r = dTr_dx / norm_Tr; float uz_r = dTr_dz / norm_Tr;            

            float nx_norm = 0.0f, nz_norm = -1.0f; 
            float cos_s = fabs(ux_s*nx_norm + uz_s*nz_norm);
            float cos_r = fabs(ux_r*nx_norm + uz_r*nz_norm);

            float a = d2Ts_dx2  + d2Tr_dx2;
            float b = d2Ts_dxdz + d2Tr_dxdz;
            float c = d2Ts_dz2  + d2Tr_dz2;

            float detH = a*c - b*b;

            float J = 1.0f / sqrtf(fabsf(detH) + EPS); 

            float R_s = max(Ts[index] / S[index], EPS);
            float R_r = max(Tr[index] / S[index], EPS);

            float G = 1.0f / sqrt(R_s * R_r);

            float theta = acosf(min(1.0f, max(-1.0f, ux_s*nx_norm + uz_s*nz_norm)));

            float R = 1.0f + 0.2f*cos(theta);

            float weights = 1.0f / (2.0f * M_PI) * sqrt(max(0.0f, cos_s) * max(0.0f, cos_r)) * G * J * R;            
            
            int mId = (i - nb) + (j - nb)*nz;
    
            atomicAdd(&data[tId], weights * model[mId]);
        }
    }    
}

__global__ void angle_domain_adjoint_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dx, float dz, float dt, float da, int nxx, int nzz, int nt, int na, int nb, int cmpId)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int i = (int)(index % nzz);
    int j = (int)(index / nzz);

    if ((i > nb) && (i < nzz-nb) && (j > nb) && (j < nxx-nb))
    {
        float T = Ts[index] + Tr[index]; 
        
        int tId = __float2int_rd(T / dt);

        int nz = (nzz - 2*nb);
        
        if (tId < nt) 
        {
            float dTs_dx = 0.5f*(Ts[i + (j+1)*nzz] - Ts[i + (j-1)*nzz]) / dx;
            float dTs_dz = 0.5f*(Ts[(i+1) + j*nzz] - Ts[(i-1) + j*nzz]) / dz;

            float d2Ts_dx2 = (Ts[i + (j+1)*nzz] - 2.0f*Ts[index] + Ts[i + (j-1)*nzz]) / (dx*dx); 
            float d2Ts_dz2 = (Ts[(i+1) + j*nzz] - 2.0f*Ts[index] + Ts[(i-1) + j*nzz]) / (dz*dz);
           
            float d2Ts_dxdz = (Ts[(i+1) + (j+1)*nzz] - Ts[(i-1) + (j+1)*nzz] - Ts[(i+1) + (j-1)*nzz] + Ts[(i-1) + (j-1)*nzz]) / (4.0f*dx*dz);

            float dTr_dx = 0.5f*(Tr[i + (j+1)*nzz] - Tr[i + (j-1)*nzz]) / dx;
            float dTr_dz = 0.5f*(Tr[(i+1) + j*nzz] - Tr[(i-1) + j*nzz]) / dz;

            float d2Tr_dx2 = (Tr[i + (j+1)*nzz] - 2.0f*Tr[index] + Tr[i + (j-1)*nzz]) / (dx*dx); 
            float d2Tr_dz2 = (Tr[(i+1) + j*nzz] - 2.0f*Tr[index] + Tr[(i-1) + j*nzz]) / (dz*dz);
           
            float d2Tr_dxdz = (Tr[(i+1) + (j+1)*nzz] - Tr[(i-1) + (j+1)*nzz] - Tr[(i+1) + (j-1)*nzz] + Tr[(i-1) + (j-1)*nzz]) / (4.0f*dx*dz);

            float norm_Ts = sqrtf(dTs_dx*dTs_dx + dTs_dz*dTs_dz) + EPS;
            float norm_Tr = sqrtf(dTr_dx*dTr_dx + dTr_dz*dTr_dz) + EPS;
            
            float ux_s = dTs_dx / norm_Ts; float uz_s = dTs_dz / norm_Ts;            
            float ux_r = dTr_dx / norm_Tr; float uz_r = dTr_dz / norm_Tr;            

            float nx_norm = 0.0f, nz_norm = -1.0f; 
            float cos_s = fabs(ux_s*nx_norm + uz_s*nz_norm);
            float cos_r = fabs(ux_r*nx_norm + uz_r*nz_norm);

            float a = d2Ts_dx2  + d2Tr_dx2;
            float b = d2Ts_dxdz + d2Tr_dxdz;
            float c = d2Ts_dz2  + d2Tr_dz2;

            float detH = a*c - b*b;

            float J = 1.0f / sqrtf(fabsf(detH) + EPS); 

            float R_s = max(Ts[index] / S[index], EPS);
            float R_r = max(Tr[index] / S[index], EPS);

            float G = 1.0f / sqrt(R_s * R_r);

            float theta = acosf(min(1.0f, max(-1.0f, ux_s*nx_norm + uz_s*nz_norm)));

            float R = 1.0f + 0.2f*cos(theta);

            float weights = 1.0f / (2.0f * M_PI) * sqrt(max(0.0f, cos_s) * max(0.0f, cos_r)) * G * J * R;            

            float reflection_angle = 0.5f*acos((dTs_dx*dTr_dx + dTs_dz*dTr_dz) / (norm_Ts * norm_Tr));

            float ang = 180.0f*reflection_angle / M_PI;
            int aId = __float2int_rd(ang / da);  

            int mId = (i - nb) + aId*nz + cmpId*na*nz;
            
            if ((aId >= 0) && (aId < na))             
                atomicAdd(&model[mId], weights * data[tId]);
        }
    }   
}

__global__ void angle_domain_forward_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dx, float dz, float dt, float da, int nxx, int nzz, int nt, int na, int nb, int cmpId)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int i = (int)(index % nzz);
    int j = (int)(index / nzz);

    if ((i > nb) && (i < nzz-nb) && (j > nb) && (j < nxx-nb))
    {
        float T = Ts[index] + Tr[index]; 
        
        int tId = __float2int_rd(T / dt);

        int nz = (nzz - 2*nb);
        
        if (tId < nt) 
        {
            float dTs_dx = 0.5f*(Ts[i + (j+1)*nzz] - Ts[i + (j-1)*nzz]) / dx;
            float dTs_dz = 0.5f*(Ts[(i+1) + j*nzz] - Ts[(i-1) + j*nzz]) / dz;

            float d2Ts_dx2 = (Ts[i + (j+1)*nzz] - 2.0f*Ts[index] + Ts[i + (j-1)*nzz]) / (dx*dx); 
            float d2Ts_dz2 = (Ts[(i+1) + j*nzz] - 2.0f*Ts[index] + Ts[(i-1) + j*nzz]) / (dz*dz);
           
            float d2Ts_dxdz = (Ts[(i+1) + (j+1)*nzz] - Ts[(i-1) + (j+1)*nzz] - Ts[(i+1) + (j-1)*nzz] + Ts[(i-1) + (j-1)*nzz]) / (4.0f*dx*dz);

            float dTr_dx = 0.5f*(Tr[i + (j+1)*nzz] - Tr[i + (j-1)*nzz]) / dx;
            float dTr_dz = 0.5f*(Tr[(i+1) + j*nzz] - Tr[(i-1) + j*nzz]) / dz;

            float d2Tr_dx2 = (Tr[i + (j+1)*nzz] - 2.0f*Tr[index] + Tr[i + (j-1)*nzz]) / (dx*dx); 
            float d2Tr_dz2 = (Tr[(i+1) + j*nzz] - 2.0f*Tr[index] + Tr[(i-1) + j*nzz]) / (dz*dz);
           
            float d2Tr_dxdz = (Tr[(i+1) + (j+1)*nzz] - Tr[(i-1) + (j+1)*nzz] - Tr[(i+1) + (j-1)*nzz] + Tr[(i-1) + (j-1)*nzz]) / (4.0f*dx*dz);

            float norm_Ts = sqrtf(dTs_dx*dTs_dx + dTs_dz*dTs_dz) + EPS;
            float norm_Tr = sqrtf(dTr_dx*dTr_dx + dTr_dz*dTr_dz) + EPS;
            
            float ux_s = dTs_dx / norm_Ts; float uz_s = dTs_dz / norm_Ts;            
            float ux_r = dTr_dx / norm_Tr; float uz_r = dTr_dz / norm_Tr;            

            float nx_norm = 0.0f, nz_norm = -1.0f; 
            float cos_s = fabs(ux_s*nx_norm + uz_s*nz_norm);
            float cos_r = fabs(ux_r*nx_norm + uz_r*nz_norm);

            float a = d2Ts_dx2  + d2Tr_dx2;
            float b = d2Ts_dxdz + d2Tr_dxdz;
            float c = d2Ts_dz2  + d2Tr_dz2;

            float detH = a*c - b*b;

            float J = 1.0f / sqrtf(fabsf(detH) + EPS); 

            float R_s = max(Ts[index] / S[index], EPS);
            float R_r = max(Tr[index] / S[index], EPS);

            float G = 1.0f / sqrt(R_s * R_r);

            float theta = acosf(min(1.0f, max(-1.0f, ux_s*nx_norm + uz_s*nz_norm)));

            float R = 1.0f + 0.2f*cos(theta);

            float weights = 1.0f / (2.0f * M_PI) * sqrt(max(0.0f, cos_s) * max(0.0f, cos_r)) * G * J * R;            

            float reflection_angle = 0.5f*acos((dTs_dx*dTr_dx + dTs_dz*dTr_dz) / (norm_Ts * norm_Tr));

            float ang = 180.0f*reflection_angle / M_PI;
            int aId = __float2int_rd(ang / da);  

            int mId = (i - nb) + aId*nz + cmpId*na*nz;
            
            if ((aId >= 0) && (aId < na))             
                atomicAdd(&data[tId], weights * model[mId]);
        }
    }   
}



