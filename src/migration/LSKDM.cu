# include "LSKDM.cuh"

void LSKDM::kirchhoff_depth_migration()
{
    prepare_convolution();

    get_src_travel_times();
    get_rec_travel_times();

    initialization();

    while (true)
    {
        compute_gradient();        
        compute_residuals();

        if (converged) break;

        compute_direction();
        compute_stepLength();

        update_model();
    }
}

void LSKDM::initialization()
{
    cal_data = new float[nt*nrx*nsy]();
    h_gradient = new float[m_samples]();
    h_direction = new float[m_samples]();
    gradient_old = new float[m_samples]();

    cudaMalloc((void**)&(d_gradient), m_samples*sizeof(float));    
    cudaMalloc((void**)&(d_direction), m_samples*sizeof(float));
    
    iteration = 0;

    cudaMemset(d_model, 0.0f, m_samples*sizeof(float));
}

void LSKDM::compute_gradient()
{
    cudaMemset(d_gradient, 0.0f, m_samples*sizeof(float));
    cudaMemcpy(d_model, h_model, m_samples*sizeof(float), cudaMemcpyHostToDevice);

    if (iteration == max_it) ++iteration;

    residuals = 0.0f;

    modeling->current_operation = domain + " LS-Migration: Computing Gradient";

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
            
            show_information();

            cudaMemcpy(d_Tr, h_Tr, nrx*volBytes, cudaMemcpyHostToDevice);
            cudaMemcpy(d_xrec, h_xrec, nrx*sizeof(float), cudaMemcpyHostToDevice);
            cudaMemcpy(d_yrec, h_yrec, nrx*sizeof(float), cudaMemcpyHostToDevice);
            
            cudaMemset(d_data, 0.0f, nt*nrx*nsy*sizeof(float));
            
            perform_forward();

            cudaMemcpy(cal_data, d_data, nt*nrx*nsy*sizeof(float), cudaMemcpyDeviceToHost);
            
            forward_convolution();
            
            for (int index = 0; index < nt*nrx*nsy; index++)
            {
                float dobs = h_data[index];
                h_data[index] = (fabsf(dobs) > EPS) ? dobs - cal_data[index] : 0.0f;
                residuals += h_data[index]*h_data[index];
            }

            adjoint_convolution();
     
            cudaMemcpy(d_data, h_data, nt*nrx*nsy*sizeof(float), cudaMemcpyHostToDevice);
            
            perform_adjoint_gradient();            
        }
    }

    if (iteration == max_it) --iteration;
    
    cudaMemcpy(h_gradient, d_gradient, m_samples*sizeof(float), cudaMemcpyDeviceToHost);
}

void LSKDM::compute_residuals()
{
    ++iteration;
    
    residuo.push_back(sqrtf(residuals));

    converged = (iteration > max_it) ? true : false;

    if (converged) std::cout << "Final residuo: "<< residuo.back() <<"\n\n";    
}

void LSKDM::compute_direction()
{   
    beta_num = 0.0f;
    beta_den = 0.0f;

    alpha_num = 0.0f;
    alpha_den = 0.0f;

    for (int index = 0; index < m_samples; index++)
    {
        beta_num += h_gradient[index]*h_gradient[index]; 
        beta_den += gradient_old[index]*gradient_old[index];
    }

    beta = beta_num / (beta_den + EPS);

    for (int index = 0; index < m_samples; index++)
    {
        h_direction[index] = beta*h_direction[index] - h_gradient[index];

        alpha_num += h_direction[index]*h_gradient[index];
    }
}

void LSKDM::compute_stepLength()
{
    cudaMemcpy(d_model, h_model, m_samples*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_direction, h_direction, m_samples*sizeof(float), cudaMemcpyHostToDevice);

    alpha = 0.0f;

    modeling->current_operation = domain + " LS-Migration: Computing Step Length";

    for (src_csIdx = 0; src_csIdx < nsx; src_csIdx++)
    {
        set_src_line_components();

        cudaMemcpy(d_Ts, h_Ts, nsy * volBytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_xsrc, h_xsrc, nsy * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_ysrc, h_ysrc, nsy * sizeof(float), cudaMemcpyHostToDevice);
        
        for (rec_csIdy = 0; rec_csIdy < nry; rec_csIdy++)
        {               
            set_rec_line_components();
            
            show_information();

            cudaMemcpy(d_Tr, h_Tr, nrx*volBytes, cudaMemcpyHostToDevice);
            cudaMemcpy(d_xrec, h_xrec, nrx*sizeof(float), cudaMemcpyHostToDevice);
            cudaMemcpy(d_yrec, h_yrec, nrx*sizeof(float), cudaMemcpyHostToDevice);
            
            cudaMemset(d_data, 0.0f, nt*nrx*nsy*sizeof(float));
            
            perform_forward_direction();

            cudaMemcpy(h_data, d_data, nt*nrx*nsy*sizeof(float), cudaMemcpyDeviceToHost);

            forward_convolution();

            for (int index = 0; index < nt*nrx*nsy; index++)
                alpha_den += h_data[index]*h_data[index];
        }
    }

    alpha = alpha_num / (alpha_den + EPS);
}

void LSKDM::update_model()
{   
    for (int index = 0; index < m_samples; index++)
    {
        h_model[index] = h_model[index] + alpha*h_direction[index];
       
        gradient_old[index] = h_gradient[index];
    }   
}

void LSKDM::show_iteration_info()
{
    if (iteration > max_it) 
        std::cout << "\n-------- Checking final residuo --------\n\n";
    else
    {    
        if (iteration == 0) 
            std::cout << "\n-------- Computing first residuo --------\n";        
        else
        {
            std::cout << "\n-------- Computing iteration " << iteration << " of " << max_it << " --------\n\n";
            
            std::cout << "Previous residuo: " << residuo.back() << "\n\n";   
        }
    }
}

void LSKDM::export_outputs()
{
    std::string path = residuo_folder + migType + "_convergence_" + std::to_string(max_it) + "_iterations.txt"; 

    std::ofstream resFile(path, std::ios::out);
    
    for (int r = 0; r < residuo.size(); r++) 
        resFile << residuo[r] << "\n";

    resFile.close();

    std::cout << "Text file \033[34m" << path << "\033[0;0m was successfully written." << std::endl;

    export_binary_float(output_path, h_model, m_samples);
}
