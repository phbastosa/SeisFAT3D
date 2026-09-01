# include "KDM.cuh"

void KDM::kirchhoff_depth_migration()
{
    prepare_convolution();

    get_src_travel_times();
    get_rec_travel_times();

    cudaMemset(d_model, 0.0f, m_samples*sizeof(float));    

    modeling->current_operation = domain + " Kirchhoff Depth Migration in " + modType + " media";

    for (src_csIdx = 0; src_csIdx < nsx; src_csIdx++)
    {
        set_src_line_components();

        cudaMemcpy(d_Ts, h_Ts, nsy * volBytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_xsrc, h_xsrc, nsy * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_ysrc, h_ysrc, nsy * sizeof(float), cudaMemcpyHostToDevice);

        for (rec_csIdy = 0; rec_csIdy < nry; rec_csIdy++)
        {               
            set_rec_line_components();
            set_cross_spread_data();
            adjoint_convolution();
     
            cudaMemcpy(d_Tr, h_Tr, nrx*volBytes, cudaMemcpyHostToDevice);
            cudaMemcpy(d_xrec, h_xrec, nrx*sizeof(float), cudaMemcpyHostToDevice);
            cudaMemcpy(d_yrec, h_yrec, nrx*sizeof(float), cudaMemcpyHostToDevice);
            cudaMemcpy(d_data, h_data, nt*nrx*nsy*sizeof(float), cudaMemcpyHostToDevice);
            
            show_information();
            perform_adjoint();            
        }
    }
}

void KDM::export_outputs()
{
    cudaMemcpy(h_model, d_model, m_samples*sizeof(float), cudaMemcpyDeviceToHost);

    export_binary_float(output_path, h_model, m_samples);
}
