# include "migration.cuh"

void Migration::set_parameters()
{
    nt = std::stoi(catch_parameter("time_samples", parameters));
    dt = std::stof(catch_parameter("time_spacing", parameters));

    max_offset = std::stof(catch_parameter("max_offset", parameters));    

    aperture_x = std::stof(catch_parameter("mig_aperture_x", parameters));
    aperture_y = std::stof(catch_parameter("mig_aperture_y", parameters));

    input_data_folder = catch_parameter("input_data_folder", parameters);
    input_data_prefix = catch_parameter("input_data_prefix", parameters);

    output_image_folder = catch_parameter("output_image_folder", parameters);
    output_table_folder = catch_parameter("output_table_folder", parameters);

    set_modeling_type();
    
    modeling->parameters = parameters;
    modeling->set_parameters();

    f_image = new float[modeling->nPoints]();
    h_image = new float[modeling->volsize]();
    h_seismic = new float[nt*modeling->max_spread]();

    cudaMalloc((void**)&(d_Tr), modeling->volsize*sizeof(float));
    cudaMalloc((void**)&(d_image), modeling->volsize*sizeof(float));
    cudaMalloc((void**)&(d_seismic), nt*modeling->max_spread*sizeof(float));

    nThreads = 256;
    nBlocks = (int)((modeling->volsize + nThreads - 1) / nThreads);
}

void Migration::read_seismic_data()
{
    std::string data_path = input_data_folder + input_data_prefix + std::to_string(modeling->geometry->sInd[modeling->srcId]+1) + ".bin";

    import_binary_float(data_path, h_seismic, nt*modeling->geometry->spread[modeling->srcId]);

    cudaMemcpy(d_seismic, h_seismic, nt*modeling->geometry->spread[modeling->srcId]*sizeof(float), cudaMemcpyHostToDevice);
}

void Migration::image_building()
{
    get_receiver_eikonal();
    run_cross_correlation();
}

void Migration::get_receiver_eikonal()
{
    for (modeling->recId = 0; modeling->recId < modeling->geometry->nrec; modeling->recId++)
    {
        set_receiver_point();

        show_information();

        modeling->time_propagation();
        
        export_receiver_eikonal();
    }
}

void Migration::set_receiver_point()
{
    modeling->sx = modeling->geometry->xrec[modeling->recId];
    modeling->sy = modeling->geometry->yrec[modeling->recId];
    modeling->sz = modeling->geometry->zrec[modeling->recId];
}

void Migration::show_information()
{
    auto clear = system("clear");
    
    std::cout << "-------------------------------------------------------------------------------\n";
    std::cout << "                                 \033[34mSeisFAT3D\033[0;0m\n";
    std::cout << "-------------------------------------------------------------------------------\n\n";

    std::cout << "Model dimensions: (z = " << (modeling->nz - 1)*modeling->dz << 
                                  ", x = " << (modeling->nx - 1)*modeling->dx <<
                                  ", y = " << (modeling->ny - 1)*modeling->dy << ") m\n\n";

    std::cout << "Running receiver " << modeling->recId + 1 << " of " << modeling->geometry->nrec << " in total\n\n";

    std::cout << "Current receiver position: (z = " << modeling->geometry->zrec[modeling->recId] << 
                                           ", x = " << modeling->geometry->xrec[modeling->recId] << 
                                           ", y = " << modeling->geometry->yrec[modeling->recId] << ") m\n\n";

    std::cout << "Kirchhoff Depth Migration: computing receiver travel time volumes\n";
}

void Migration::export_receiver_eikonal()
{
    cudaMemcpy(modeling->T, modeling->d_T, modeling->volsize*sizeof(float), cudaMemcpyDeviceToHost);

    export_binary_float(output_table_folder + "eikonal_receiver_" + std::to_string(modeling->recId+1) + ".bin", modeling->T, modeling->nPoints);    
}

void Migration::run_cross_correlation()
{
    cudaMemset(d_image, 0.0f, modeling->volsize*sizeof(float));

    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nrel; modeling->srcId++)
    {
        read_seismic_data();

        modeling->set_shot_point();
        modeling->show_information();
        modeling->time_propagation();

        std::cout << "\nKirchhoff depth migration: computing image matrix\n";

        int spread = 0;
        
        for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
        {
            float rx = modeling->geometry->xrec[modeling->recId];
            float ry = modeling->geometry->yrec[modeling->recId];

            float offset = sqrtf(powf(modeling->sx - rx, 2.0f) + powf(modeling->sy - ry, 2.0f));

            float cmp_x = modeling->sx + 0.5f*(rx - modeling->sx);
            float cmp_y = modeling->sy + 0.5f*(ry - modeling->sy);

            if (offset < max_offset)
            {
                import_binary_float(output_table_folder + "eikonal_receiver_" + std::to_string(modeling->recId+1) + ".bin", modeling->T, modeling->volsize);
            
                cudaMemcpy(d_Tr, modeling->T, modeling->volsize*sizeof(float), cudaMemcpyHostToDevice);

                cross_correlation<<<nBlocks, nThreads>>>(modeling->d_T, d_Tr, d_image, d_seismic, aperture_x, aperture_y, cmp_x, cmp_y, spread, modeling->nxx, 
                                                         modeling->nyy, modeling->nzz, modeling->nb, modeling->dx, modeling->dy, modeling->dz, nt, dt);
            }

            ++spread;
        }
    }
}

void Migration::export_outputs()
{
    cudaMemcpy(h_image, d_image, modeling->volsize*sizeof(float), cudaMemcpyDeviceToHost);
    modeling->reduce_boundary(h_image, f_image);
    export_binary_float(output_image_folder + "kirchhoff_result_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + "x" + std::to_string(modeling->ny) + ".bin", f_image, modeling->nPoints);
}

__global__ void cross_correlation(float * Ts, float * Tr, float * image, float * seismic, float aperture_x, float aperture_y, float cmp_x, 
                                  float cmp_y, int spread, int nxx, int nyy, int nzz, int nb, float dx, float dy, float dz, int nt, float dt)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         
    int j = (int) (index - k*nxx*nzz) / nzz;   
    int i = (int) (index - j*nzz - k*nxx*nzz); 

    if ((i > nb) && (i < nzz-nb) && (j > nb) && (j < nxx-nb) && (k > nb) && (k < nyy-nb))
    {
        float sigma_x = tanf(aperture_x * M_PI / 180.0f)*(i-nb)*dz;        
        float sigma_y = tanf(aperture_y * M_PI / 180.0f)*(i-nb)*dz;        

        float par_x = powf(((j-nb)*dx - cmp_x) / (sigma_x + 1e-6f), 2.0f);
        float par_y = powf(((k-nb)*dy - cmp_y) / (sigma_y + 1e-6f), 2.0f);

        float value = expf(-0.5f*(par_x + par_y));

        float T = Ts[index] + Tr[index]; 
    
        int tId = (int)(T / dt);

        if (tId < nt) image[index] += value * seismic[tId + spread*nt];
    }
}    

