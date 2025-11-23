# include "LSKDM.cuh"

void LSKDM::kirchhoff_depth_migration()
{
    set_src_travel_times();
    set_rec_travel_times();
    prepare_convolution();

    // set_src_domain();

    // initialization();

    // while (true)
    // {
    //     compute_gradient();        
    //     compute_residuals();

    //     if (converged) break;

    //     compute_direction();
    //     compute_stepLength();

    //     update_model();
    // }
}

void LSKDM::initialization()
{
    // h_gradient = new float[m_samples]();
    // h_direction = new float[m_samples]();
    // gradient_old = new float[m_samples]();

    // cudaMalloc((void**)&(d_gradient), m_samples*sizeof(float));    
    // cudaMalloc((void**)&(d_direction), m_samples*sizeof(float));

    // iteration = 0;

    // cudaMemset(d_model, 0.0f, m_samples*sizeof(float));

    // for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nsrc; modeling->srcId++)
    // {     
    //     std::string data_path = input_data_folder + input_data_prefix + std::to_string(modeling->srcId+1) + ".bin";
    //     import_binary_float(data_path, seismic, nt*modeling->max_spread);

    //     import_binary_float(tables_folder + "eikonal_src_" + std::to_string(modeling->srcId+1) + ".bin", h_Ts, modeling->matsize);
    //     cudaMemcpy(d_Ts, h_Ts, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

    //     set_current_src();

    //     current_operation = domain + " LS-Migration: Computing Initial Model";

    //     show_information();
    //     show_iteration_info();

    //     int spreadId = 0;
        
    //     for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
    //     {
    //         cmpId = spreadId + 2.0f*(ds/dr)*modeling->srcId;                              
            
    //         import_binary_float(tables_folder + "eikonal_rec_" + std::to_string(modeling->recId+1) + ".bin", h_Tr, modeling->matsize);
    //         cudaMemcpy(d_Tr, h_Tr, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

    //         for (int tId = 0; tId < nt; tId++)
    //             h_data[tId] = seismic[tId + spreadId*nt];
            
    //         adjoint_convolution();
            
    //         cudaMemcpy(d_data, h_data, nt * sizeof(float), cudaMemcpyHostToDevice);

    //         perform_adjoint();

    //         ++spreadId;
    //     }
    // }
    
    // cudaMemcpy(h_model, d_model, m_samples*sizeof(float), cudaMemcpyDeviceToHost);
}

void LSKDM::compute_gradient()
{
    // cudaMemset(d_gradient, 0.0f, m_samples*sizeof(float));
    // cudaMemcpy(d_model, h_model, m_samples*sizeof(float), cudaMemcpyHostToDevice);

    // if (iteration == max_it) ++iteration;

    // residuals = 0.0f;

    // for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nsrc; modeling->srcId++)
    // {     
    //     std::string data_path = input_data_folder + input_data_prefix + std::to_string(modeling->srcId+1) + ".bin";
    //     import_binary_float(data_path, seismic, nt*modeling->max_spread);

    //     import_binary_float(tables_folder + "eikonal_src_" + std::to_string(modeling->srcId+1) + ".bin", h_Ts, modeling->matsize);
    //     cudaMemcpy(d_Ts, h_Ts, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

    //     set_current_src();

    //     current_operation = domain + " LS-Migration: Computing Gradient";

    //     show_information();
    //     show_iteration_info();

    //     int spreadId = 0;
        
    //     for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
    //     {
    //         cmpId = spreadId + 2.0f*(ds/dr)*modeling->srcId;                              

    //         import_binary_float(tables_folder + "eikonal_rec_" + std::to_string(modeling->recId+1) + ".bin", h_Tr, modeling->matsize);
    //         cudaMemcpy(d_Tr, h_Tr, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

    //         cudaMemset(d_data, 0.0f, nt * sizeof(float));

    //         perform_forward();

    //         cudaMemcpy(h_data, d_data, nt * sizeof(float), cudaMemcpyDeviceToHost);

    //         forward_convolution();

    //         for (int tId = 0; tId < nt; tId++)
    //         {                
    //             h_data[tId] -= seismic[tId + spreadId*nt];
                
    //             residuals += h_data[tId]*h_data[tId];
    //         }
            
    //         adjoint_convolution();
            
    //         cudaMemcpy(d_data, h_data, nt * sizeof(float), cudaMemcpyHostToDevice);
            
    //         perform_adjoint_gradient();

    //         ++spreadId;
    //     }
    // }

    // if (iteration == max_it) --iteration;
    
    // cudaMemcpy(h_gradient, d_gradient, m_samples*sizeof(float), cudaMemcpyDeviceToHost);

    // float grad_norm = 0.0f;
    // for (int index = 0; index < m_samples; index++)
    //     grad_norm += h_gradient[index]*h_gradient[index];

    // grad_norm = sqrtf(grad_norm);

    // for (int index = 0; index < m_samples; index++)
    //     h_gradient[index] /= grad_norm;
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
    // float beta_num = 0.0f;
    // float beta_den = 0.0f;

    // for (int index = 0; index < m_samples; index++)
    // {
    //     beta_num += (h_gradient[index] - gradient_old[index])*h_gradient[index]; 
    //     beta_den += gradient_old[index]*gradient_old[index];
    // }

    // beta = beta_num / (beta_den + EPS);

    // for (int index = 0; index < m_samples; index++)
    //     h_direction[index] = beta*h_direction[index] - h_gradient[index];
}

void LSKDM::compute_stepLength()
{
    // cudaMemcpy(d_model, h_model, m_samples*sizeof(float), cudaMemcpyHostToDevice);
    // cudaMemcpy(d_direction, h_direction, m_samples*sizeof(float), cudaMemcpyHostToDevice);

    // alpha = 0.0f;

    // float alpha_num = 0.0f;
    // float alpha_den = 0.0f;

    // for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nsrc; modeling->srcId++)
    // {     
    //     std::string data_path = input_data_folder + input_data_prefix + std::to_string(modeling->srcId+1) + ".bin";
    //     import_binary_float(data_path, seismic, nt*modeling->max_spread);

    //     import_binary_float(tables_folder + "eikonal_src_" + std::to_string(modeling->srcId+1) + ".bin", h_Ts, modeling->matsize);
    //     cudaMemcpy(d_Ts, h_Ts, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

    //     set_current_src();

    //     current_operation = domain + " LS-Migration: Computing Step Length";

    //     show_information();
    //     show_iteration_info();
        
    //     int spreadId = 0;
        
    //     for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
    //     {
    //         cmpId = spreadId + 2.0f*(ds/dr)*modeling->srcId;                              

    //         import_binary_float(tables_folder + "eikonal_rec_" + std::to_string(modeling->recId+1) + ".bin", h_Tr, modeling->matsize);
    //         cudaMemcpy(d_Tr, h_Tr, modeling->matsize*sizeof(float), cudaMemcpyHostToDevice);

    //         cudaMemset(d_data, 0.0f, nt * sizeof(float));

    //         perform_forward();

    //         cudaMemcpy(h_data, d_data, nt * sizeof(float), cudaMemcpyDeviceToHost);

    //         forward_convolution();

    //         for (int tId = 0; tId < nt; tId++)
    //             seismic[tId + spreadId*nt] = seismic[tId + spreadId*nt] - h_data[tId];
            
    //         cudaMemset(d_data, 0.0f, nt * sizeof(float));

    //         perform_forward_direction();
            
    //         cudaMemcpy(h_data, d_data, nt * sizeof(float), cudaMemcpyDeviceToHost);

    //         forward_convolution();
            
    //         for (int tId = 0; tId < nt; tId++)
    //         {
    //             alpha_num += seismic[tId + spreadId*nt]*h_data[tId];
    //             alpha_den += h_data[tId]*h_data[tId];
    //         }    
            
    //         ++spreadId;
    //     }
    // }
    
    // alpha = alpha_num / (alpha_den + EPS);
}

void LSKDM::update_model()
{
    // for (int index = 0; index < m_samples; index++)
    // {
    //     h_model[index] = h_model[index] + alpha*h_direction[index];
        
    //     gradient_old[index] = h_gradient[index];
    // }   
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
