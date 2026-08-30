# include "KDM.cuh"

void KDM::kirchhoff_depth_migration()
{
    prepare_convolution();

    get_src_travel_times();
    get_rec_travel_times();

    cudaMemset(d_model, 0.0f, m_samples*sizeof(float));    

    modeling->current_operation = domain + " Kirchhoff Depth Migration in " + modType + " media";

    auto set_rec_line_travel_times = [&](float * buffer, int rec_csIdy) 
    {
        for (int rec_csIdx = 0; rec_csIdx < nrx; rec_csIdx++) 
        {
            int recId = rec_csIdy + rec_csIdx*nry;

            float * target = buffer + rec_csIdx*modeling->volsize;
            
            std::string path = tables_folder + "eikonal_rec_" + std::to_string(recId+1) + ".bin";
            
            import_binary_float(path, target, modeling->volsize);
        }
    };

    auto set_cross_spread_data = [&](float * buffer, int rec_csIdy) 
    {
        #pragma omp parallel for collapse(2) schedule(static)
        for (int src_csIdy = 0; src_csIdy < nsy; src_csIdy++)
        {
            for (int rec_csIdx = 0; rec_csIdx < nrx; rec_csIdx++)
            {            
                size_t d_rec_offset = rec_csIdx*nt;
                size_t d_src_offset = src_csIdy*nt*nrx;
                size_t h_rec_offset = (rec_csIdy + rec_csIdx*nry)*nt;
                size_t h_src_offset = src_csIdy*nt*modeling->geometry->nrec;

                size_t hst_base = h_rec_offset + h_src_offset; 
                size_t dvc_base = d_rec_offset + d_src_offset;

                std::memcpy(&buffer[dvc_base], &seismic[hst_base], nt*sizeof(float));
            }
        }
    };

    auto perform_adjoint_convolution = [&](float * buffer)
    {
        int num_traces = nrx * nsy; 

        for (int csId = 0; csId < num_traces; csId++)
        {
            for (int tId = 0; tId < nfft; tId++)
            {
                if (tId < nt)
                {
                    time_trace[tId] = (double)buffer[tId + csId*nt];
                    buffer[tId + csId*nt] = 0.0f;
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

            for (int tId = nw/2; tId < nt; tId++)
                buffer[tId + csId*nt] = (float)(time_trace[tId - nw/2]);
        }
    };

    auto set_rec_line_geometry = [&](float * h_xrec, float * h_yrec, int rec_csIdy) 
    {
        for (int rec_csIdx = 0; rec_csIdx < nrx; rec_csIdx++) 
        {
            int recId = rec_csIdy + rec_csIdx*nry;

            h_xrec[rec_csIdx] = modeling->geometry->xrec[recId];
            h_yrec[rec_csIdx] = modeling->geometry->yrec[recId];
        }
    };

    for (src_csIdx = 0; src_csIdx < nsx; src_csIdx++)
    {
        set_src_line_components();

        set_cross_spread_data(h_data[0], 0);
        set_rec_line_travel_times(h_Tr[0], 0);
        perform_adjoint_convolution(h_data[0]);
        set_rec_line_geometry(h_xrec[0], h_yrec[0], 0);

        for (rec_csIdy = 0; rec_csIdy < nry; rec_csIdy++)
        {               
            curr = rec_csIdy % 2;

            int next = (rec_csIdy + 1) % 2;
            int next_recId = rec_csIdy + 1;

            zpos = format1Decimal(modeling->geometry->zrec[rec_csIdy]);
            ypos = format1Decimal(modeling->geometry->xrec[rec_csIdy]);
            xpos = format1Decimal(modeling->geometry->yrec[rec_csIdy +       0*nry]) + " - " +
                   format1Decimal(modeling->geometry->yrec[rec_csIdy + (nrx-1)*nry]);

            cudaMemcpyAsync(d_data[curr], h_data[curr], nt*nrx*nsy*sizeof(float), cudaMemcpyHostToDevice, stream_cpy);
            cudaMemcpyAsync(d_Tr[curr], h_Tr[curr], nrx*volBytes, cudaMemcpyHostToDevice, stream_cpy);
            cudaMemcpyAsync(d_xrec[curr], h_xrec[curr], nrx*sizeof(float), cudaMemcpyHostToDevice, stream_cpy);
            cudaMemcpyAsync(d_yrec[curr], h_yrec[curr], nrx*sizeof(float), cudaMemcpyHostToDevice, stream_cpy);
            cudaEventRecord(cpy_done[curr], stream_cpy);            

            cudaStreamWaitEvent(stream_krn, cpy_done[curr], 0);        
            
            show_information();
            perform_adjoint();            
            
            cudaEventRecord(krn_done[curr], stream_krn);

            std::future<void> worker;

            if (rec_csIdy + 1 < nry) 
            {
                if (rec_csIdy >= 1) cudaEventSynchronize(krn_done[next]);
                                
                worker = std::async(std::launch::async, [&, next, next_recId]() 
                {
                    set_cross_spread_data(h_data[next], next_recId);
                    set_rec_line_travel_times(h_Tr[next], next_recId);
                    perform_adjoint_convolution(h_data[next]);
                    set_rec_line_geometry(h_xrec[next], h_yrec[next], next_recId);
                });
            }

            if (worker.valid()) worker.wait();    
        }
    }
}

void KDM::export_outputs()
{
    cudaMemcpy(h_model, d_model, m_samples*sizeof(float), cudaMemcpyDeviceToHost);

    export_binary_float(output_path, h_model, m_samples);
}
