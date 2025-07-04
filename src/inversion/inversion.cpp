# include "inversion.hpp"

void Inversion::set_parameters()
{
    max_iteration = std::stoi(catch_parameter("max_iteration", parameters));

    tk_order = std::stoi(catch_parameter("tk_order", parameters));
    tk_param = std::stof(catch_parameter("tk_param", parameters));

    obs_data_folder = catch_parameter("obs_data_folder", parameters);
    obs_data_prefix = catch_parameter("obs_data_prefix", parameters);

    smoother_stdv = std::stoi(catch_parameter("gaussian_filter_stdv", parameters));
    smoother_samples = std::stoi(catch_parameter("gaussian_filter_samples", parameters));
    
    convergence_map_folder = catch_parameter("convergence_folder", parameters);
    estimated_model_folder = catch_parameter("inversion_output_folder", parameters);

    smooth_model_per_iteration = str2bool(catch_parameter("smooth_per_iteration", parameters));

    set_modeling_type();

    n_model = modeling->nPoints;

    x = new float[n_model]();
}

void Inversion::import_obsData()
{
    n_data = modeling->max_spread * modeling->geometry->nrel;

    dcal = new float[n_data]();
    dobs = new float[n_data]();

    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nrel; modeling->srcId++)
    {
        float * data = new float[modeling->max_spread]();

        std::string path = obs_data_folder + obs_data_prefix + std::to_string(modeling->geometry->sInd[modeling->srcId]+1) + ".bin";

        import_binary_float(path, data, modeling->max_spread);

        int skipped = modeling->srcId * modeling->max_spread;    
        
        for (int i = 0; i < modeling->max_spread; i++) 
            dobs[i + skipped] = data[i];

        delete[] data;
    }
}

void Inversion::forward_modeling()
{
    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nrel; modeling->srcId++)
    {
        modeling->set_shot_point();
        
        show_information();

        modeling->time_propagation();
        
        concatenate_data();
        
        if (iteration != max_iteration)
            gradient_ray_tracing();
    }
}

void Inversion::show_information()
{
    modeling->show_information();    
    
    std::cout << "\nInversion type: " << inversion_method << "\n\n";

    if (iteration == max_iteration) 
        std::cout << "-------- Checking final residuo --------\n\n";
    else
    {    
        std::cout << "-------- Computing iteration " << iteration + 1 << " of " << max_iteration << " --------\n\n";

        if (iteration > 0) std::cout << "Previous residuo: " << residuo.back() << "\n\n";   
    }
}

void Inversion::concatenate_data()
{
    modeling->compute_seismogram();

    int skipped = modeling->srcId * modeling->max_spread;

    for (int i = 0; i < modeling->max_spread; i++) 
        dcal[i + skipped] = modeling->seismogram[i];    
}

void Inversion::gradient_ray_tracing()
{
    int sIdx = (int)((modeling->sx + 0.5f*modeling->dx) / modeling->dx);
    int sIdy = (int)((modeling->sy + 0.5f*modeling->dy) / modeling->dy);
    int sIdz = (int)((modeling->sz + 0.5f*modeling->dz) / modeling->dz);

    int sId = sIdz + sIdx*modeling->nz + sIdy*modeling->nx*modeling->nz; 

    float rayStep = 0.2f*modeling->dz;

    std::vector<int> ray_index; 

    for (int ray_id = modeling->geometry->iRec[modeling->srcId]; ray_id < modeling->geometry->fRec[modeling->srcId]; ray_id++)
    {
        float xi = modeling->geometry->xrec[ray_id];        
        float yi = modeling->geometry->yrec[ray_id];        
        float zi = modeling->geometry->zrec[ray_id];

        if ((modeling->sz == zi) && (modeling->sx == xi) && (modeling->sy == zi)) continue;        

        while (true)
        {
            int k = (int)((yi + 0.5f*modeling->dy) / modeling->dy) + modeling->nb; 
            int j = (int)((xi + 0.5f*modeling->dx) / modeling->dx) + modeling->nb; 
            int i = (int)((zi + 0.5f*modeling->dz) / modeling->dz) + modeling->nb; 

            float dTz = 0.5f*(modeling->T[(i+1) + j*modeling->nzz + k*modeling->nxx*modeling->nzz] - modeling->T[(i-1) + j*modeling->nzz + k*modeling->nxx*modeling->nzz]) / modeling->dz;    
            float dTx = 0.5f*(modeling->T[i + (j+1)*modeling->nzz + k*modeling->nxx*modeling->nzz] - modeling->T[i + (j-1)*modeling->nzz + k*modeling->nxx*modeling->nzz]) / modeling->dx;    
            float dTy = 0.5f*(modeling->T[i + j*modeling->nzz + (k+1)*modeling->nxx*modeling->nzz] - modeling->T[i + j*modeling->nzz + (k-1)*modeling->nxx*modeling->nzz]) / modeling->dy;    

            float norm = sqrtf(dTx*dTx + dTy*dTy + dTz*dTz);

            xi -= rayStep*dTx / norm;   
            yi -= rayStep*dTy / norm;   
            zi -= rayStep*dTz / norm;    

            int km = (int)((yi + 0.5f*modeling->dy) / modeling->dy); 
            int jm = (int)((xi + 0.5f*modeling->dx) / modeling->dx); 
            int im = (int)((zi + 0.5f*modeling->dz) / modeling->dz); 
            
            int index = im + jm*modeling->nz + km*modeling->nx*modeling->nz;

            ray_index.push_back(index);

            if (ray_index.back() == sId) break;
        }
   
        float final_distance = sqrtf(powf(zi - modeling->sz, 2.0f) + 
                                     powf(xi - modeling->sx, 2.0f) + 
                                     powf(yi - modeling->sy, 2.0f));

        std::sort(ray_index.begin(), ray_index.end());

        int current_voxel_index = ray_index[0];
        float distance_per_voxel = rayStep;

        for (int index = 0; index < ray_index.size(); index++)
        {
            if (ray_index[index] == current_voxel_index)
            {
                distance_per_voxel += rayStep;
            }
            else
            {
                vG.push_back(distance_per_voxel);
                jG.push_back(current_voxel_index);
                iG.push_back(ray_id + modeling->srcId * modeling->max_spread);

                if (current_voxel_index == sId) vG.back() = final_distance;

                distance_per_voxel = rayStep;
                current_voxel_index = ray_index[index];    
            }
        }

        if (current_voxel_index == sId)
        {
            vG.push_back(final_distance);
            jG.push_back(current_voxel_index);
            iG.push_back(ray_id + modeling->srcId * modeling->max_spread);
        }
        else 
        {
            vG.push_back(distance_per_voxel);
            jG.push_back(current_voxel_index);
            iG.push_back(ray_id + modeling->srcId * modeling->max_spread);
        }

        std::vector<int>().swap(ray_index);
    }
}

void Inversion::check_convergence()
{
    set_objective_function();

    if ((iteration >= max_iteration))
    {
        std::cout << "Final residuo: "<< residuo.back() <<"\n";
        converged = true;
    }
    else
    {
        iteration += 1;
        converged = false;
    }
}

void Inversion::optimization()
{
    set_sensitivity_matrix();
    
    set_regularization();
    
    solve_linear_system_lscg();

    delete[] iA;
    delete[] jA;
    delete[] vA;
    delete[] B;
}

void Inversion::solve_linear_system_lscg()
{
    float a, b, qTq, rTr, rd;
    int cg_max_iteration = 10;

    float * s = new float[N]();
    float * q = new float[N]();
    float * r = new float[M]();
    float * p = new float[M]();

    for (int i = 0; i < N; i++) 
        s[i] = B[i]; 

    for (int i = 0; i < NNZ; i++) 
        r[jA[i]] += vA[i] * s[iA[i]];        

    for (int i = 0; i < M; i++) 
    {
        p[i] = r[i]; 
        x[i] = 0.0f;
    }

    for (int i = 0; i < NNZ; i++) 
        q[iA[i]] += vA[i] * p[jA[i]];        

    for (int i = 0; i < cg_max_iteration; i++)
    {
        qTq = 0.0f;
        for (int k = 0; k < N; k++)           
            qTq += q[k] * q[k];               

        rTr = 0.0f;
        for (int k = 0; k < M; k++)           
            rTr += r[k] * r[k];                

        a = rTr / qTq;                                            

        for (int k = 0; k < M; k++)           
            x[k] += a * p[k];                 

        for (int k = 0; k < N; k++)             
            s[k] -= a * q[k];                  

        rd = 0.0f;
        for (int k = 0; k < M; k++)            
            rd += r[k] * r[k];                

        for (int k = 0; k < M; k++)            
            r[k] = 0.0f;                      
        
        for (int k = 0; k < NNZ; k++)          
            r[jA[k]] += vA[k] * s[iA[k]];         

        rTr = 0.0f;                
        for (int k = 0; k < M; k++)           
            rTr += r[k] * r[k];               

        b = rTr / rd;                         

        for (int k = 0; k < M; k++)          
            p[k] = r[k] + b * p[k];            

        for (int k = 0; k < N; k++) 
            q[k] = 0.0f;                      

        for (int k = 0; k < NNZ; k++) 
            q[iA[k]] += vA[k] * p[jA[k]];      
    }

    delete[] s;
    delete[] q;
    delete[] r;
    delete[] p;
}

void Inversion::model_update()
{
    if (smooth_model_per_iteration)
    {
        int aux_nx = modeling->nx + 2*smoother_samples;
        int aux_ny = modeling->ny + 2*smoother_samples;
        int aux_nz = modeling->nz + 2*smoother_samples;

        int aux_nPoints = aux_nx*aux_ny*aux_nz;

        float * dm_aux = new float[aux_nPoints]();
        float * dm_smooth = new float[aux_nPoints]();

        # pragma omp parallel for
        for (int index = 0; index < modeling->nPoints; index++)
        {
            int k = (int) (index / (modeling->nx*modeling->nz));         
            int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;   
            int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz); 

            int ind_filt = (i + smoother_samples) + (j + smoother_samples)*aux_nz + (k + smoother_samples)*aux_nx*aux_nz;

            dm_aux[ind_filt] = x[i + j*modeling->nz + k*modeling->nx*modeling->nz];
        }

        smooth_volume(dm_aux, dm_smooth, aux_nx, aux_ny, aux_nz);

        # pragma omp parallel for    
        for (int index = 0; index < modeling->nPoints; index++)
        {
            int k = (int) (index / (modeling->nx*modeling->nz));         
            int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;   
            int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz); 

            int ind_filt = (i + smoother_samples) + (j + smoother_samples)*aux_nz + (k + smoother_samples)*aux_nx*aux_nz;

            x[i + j*modeling->nz + k*modeling->nx*modeling->nz] = dm_smooth[ind_filt];
        }
    
        delete[] dm_aux;
        delete[] dm_smooth;
    }   
    
    for (int index = 0; index < modeling->nPoints; index++)
    {
        int k = (int) (index / (modeling->nx*modeling->nz));         
        int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;   
        int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz); 

        int indb = (i + modeling->nb) + (j + modeling->nb)*modeling->nzz + (k + modeling->nb)*modeling->nxx*modeling->nzz;

        modeling->S[indb] += x[index];
    }

    modeling->copy_slowness_to_device();
}

void Inversion::smooth_volume(float * input, float * output, int nx, int ny, int nz)
{
    int nPoints = nx * ny * nz;
    int nKernel = smoother_samples * smoother_samples * smoother_samples;

    float * kernel = new float[nKernel]();

    # pragma omp parallel for
    for (int i = 0; i < nPoints; i++) 
        output[i] = input[i];

    int mid = (int)(smoother_samples / 2); 

    kernel[mid + mid*smoother_samples + mid*smoother_samples*smoother_samples] = 1.0f;

    if (smoother_stdv != 0.0f)
    {
        float sum = 0.0f;

        for (int y = -mid; y <= mid; y++)
        {
            for (int x = -mid; x <= mid; x++)
            {
                for (int z = -mid; z <= mid; z++)
                {          
                    int index = (z + mid) + (x + mid)*smoother_samples + (y + mid)*smoother_samples*smoother_samples; 
                    
                    float r = sqrtf(x*x + y*y + z*z);

                    kernel[index] = 1.0f / (M_PI*smoother_stdv) * expf(-((r*r)/(2.0f*smoother_stdv*smoother_stdv)));
        
                    sum += kernel[index]; 
                }
            }
        }

        for (int i = 0; i < nKernel; i++) 
            kernel[i] /= sum;
    }
        
    for (int k = mid; k < ny - mid; k++)
    {   
        for (int j = mid; j < nx - mid; j++)
        {
            for (int i = mid; i < nz - mid; i++)
            {       
                float accum = 0.0f;
                
                for (int yk = 0; yk < smoother_samples; yk++)
                {      
                    for (int xk = 0; xk < smoother_samples; xk++)
                    {      
                        for (int zk = 0; zk < smoother_samples; zk++)
                        {   
                            int index = zk + xk*smoother_samples + yk*smoother_samples*smoother_samples;   
                            int partial = (i - mid + zk) + (j - mid + xk)*nz + (k - mid + yk)*nx*nz; 

                            accum += input[partial] * kernel[index];
                        }        
                    }
                }
                
                output[i + j*nz + k*nx*nz] = accum;
            }
        }   
    }

    delete[] kernel;
}

void Inversion::export_results()
{    
    std::string estimated_model_path = estimated_model_folder + inversion_name + "_final_model_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + "x" + std::to_string(modeling->ny) + ".bin";
    std::string convergence_map_path = convergence_map_folder + inversion_name + "_convergence_" + std::to_string(iteration) + "_iterations.txt"; 

    float * Vp = new float[modeling->nPoints]();
    modeling->reduce_boundary(modeling->S, Vp);

    # pragma omp parallel for
    for (int index = 0; index < modeling->nPoints; index++)
        Vp[index] = 1.0f / Vp[index];

    export_binary_float(estimated_model_path, Vp, modeling->nPoints);

    std::ofstream resFile(convergence_map_path, std::ios::out);
    
    for (int r = 0; r < residuo.size(); r++) 
        resFile << residuo[r] << "\n";

    resFile.close();

    std::cout << "Text file \033[34m" << convergence_map_path << "\033[0;0m was successfully written." << std::endl;
}
