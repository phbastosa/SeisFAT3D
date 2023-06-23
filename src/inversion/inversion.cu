# include "inversion.cuh"

void Inversion::set_parameters()
{
    Modeling::file = file;
    Modeling::set_parameters();

    iteration = 0;

    n_model = nPoints;
    n_data = total_nodes * total_shots;

    tolerance = std::stof(catch_parameter("tolerance", file));
    max_iteration = std::stof(catch_parameter("max_iteration", file));

    obs_data_folder = catch_parameter("obs_data_folder", file);
    obs_data_prefix = catch_parameter("obs_data_prefix", file);

    dobs = new float[n_data]();    
    dcal = new float[n_data]();    

    source = new float[n_model]();
    adjoint = new float[n_model]();
    gradient = new float[n_model]();
}

void Inversion::import_obs_data()
{    
    int ptr = 0; 

    float * data = new float[total_nodes]();

    for (int shot = 0; shot < total_shots; shot++)
    {
        import_binary_float(obs_data_folder + obs_data_prefix + std::to_string(shot+1) + ".bin", data, total_nodes);

        for (int node = ptr; node < ptr + total_nodes; node++) 
            dobs[node] = data[node - ptr];

        ptr += total_nodes;        
    }

    delete[] data;

    set_runtime();
}

void Inversion::forward_modeling()
{
    set_slowness();

    for (int shot = 0; shot < total_shots; shot++)
    {
        shot_id = shot;

        info_message();

        initial_setup();

        forward_solver();

        build_outputs();

        fill_calculated_data();

        adjoint_state_solver();



    }

    free_space(); 
}

void Inversion::fill_calculated_data()
{
    for (int node = 0; node < total_nodes; node++)
    {
        int index = node + shot_id * total_nodes;

        dcal[index] = receiver_output[node]; 
    }
}

void Inversion::adjoint_state_solver()
{
    float cell_volume = dx * dy * dz;

    for (int index = 0; index < n_model; index++) 
    {
        source[index] = 0.0f;    
        adjoint[index] = 1e6f;

        int k = (int) (index / (nx*nz));        
        int j = (int) (index - k*nx*nz) / nz;    
        int i = (int) (index - j*nz - k*nx*nz);  

        if ((i == 0) || (i == nz-1) || 
            (j == 0) || (j == nx-1) || 
            (k == 0) || (k == ny-1))  
        {    
            adjoint[index] = 0.0f;        
        }
    }

    for (int node = 0; node < total_nodes; node++)
    {
        int node_id = node + shot_id * total_nodes;

        float tmp = dcal[node_id] + t0 - dobs[node_id]; 

        int i = (int) floorf(geometry->nodes.z[node] / dz);
        int j = (int) floorf(geometry->nodes.x[node] / dx);
        int k = (int) floorf(geometry->nodes.y[node] / dy);

        int index = i + j*nz + k*nz*nx;

        source[index] += tmp / cell_volume; 
        source[index + 1] += tmp / cell_volume; 
        source[index + nz] += tmp / cell_volume; 
        source[index + 1 + nz] += tmp / cell_volume; 
        source[index + nx * nz] += tmp / cell_volume; 
        source[index + 1 + nx * nz] += tmp / cell_volume; 
        source[index + nz + nx * nz] += tmp / cell_volume; 
        source[index + 1 + nz + nx * nz] += tmp / cell_volume; 
    }




    for (int sweep = 0; sweep < nSweeps; sweep++)
	{ 
		int start = (sweep == 3 || sweep == 5 || sweep == 6 || sweep == 7) ? totalLevels : meshDim;
		int end = (start == meshDim) ? totalLevels + 1 : meshDim - 1;
		int incr = (start == meshDim) ? true : false;

		int xSweepOff = (sweep == 3 || sweep == 4) ? nxx : 0;
		int ySweepOff = (sweep == 2 || sweep == 5) ? nyy : 0;
		int zSweepOff = (sweep == 1 || sweep == 6) ? nzz : 0;
		
		for (int level = start; level != end; level = (incr) ? level + 1 : level - 1)
		{			
			int xs = max(1, level - (nyy + nzz));	
			int ys = max(1, level - (nxx + nzz));

			int xe = min(nxx, level - (meshDim - 1));
			int ye = min(nyy, level - (meshDim - 1));	
		
			int xr = xe - xs + 1;
			int yr = ye - ys + 1;

			int nThreads = xr * yr;
				
			dim3 bs(16, 16, 1);

			if (nThreads < threadsPerBlock) { bs.x = xr; bs.y = yr; } 

			dim3 gs(iDivUp(xr, bs.x), iDivUp(yr , bs.y), 1);
			
            int sgni = sweep + 0*nSweeps;
            int sgnj = sweep + 1*nSweeps;
            int sgnk = sweep + 2*nSweeps;


            // adjoint_state_kernel();
            // cudaDeviceSyncronize();
		}
	}









    // for (int index = 0; index < n_model; index++) 
    // {
    //     float slowness = 1.0f / V[index];

    //     gradient[index] += source[index];
    // }
}

void Inversion::compute_residuals()
{
    float r = 0.0f;
    for (int i = 0; i < n_data; i++)
        r += powf(dobs[i] - dcal[i], 2.0f);

    residuo.push_back(sqrtf(r));
}

void Inversion::check_convergence()
{
    if ((iteration >= max_iteration) || (residuo.back() <= tolerance))
    {
        std::cout<<"\nFinal residuo: "<<residuo.back()<<std::endl;
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


}

void Inversion::model_update()
{



}



void Inversion::export_results()
{
    get_runtime();


    export_binary_float("gradient.bin", gradient, n_model);
}