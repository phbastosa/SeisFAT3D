# include "tomography.cuh"

void Tomography::set_parameters()
{
    modeling = new Eikonal();

    modeling->file = file;

    modeling->set_parameters();

    iteration = 0;

    n_model = modeling->nPoints;
    n_data = modeling->total_nodes * modeling->total_shots;

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

void Tomography::import_obs_data()
{    
    int ptr = 0; 

    float * data = new float[modeling->total_nodes]();

    for (int shot = 0; shot < modeling->total_shots; shot++)
    {
        import_binary_float(obs_data_folder + obs_data_prefix + std::to_string(shot+1) + ".bin", data, modeling->total_nodes);

        for (int node = ptr; node < ptr + modeling->total_nodes; node++) 
            dobs[node] = data[node - ptr];

        ptr += modeling->total_nodes;        
    }

    delete[] data;

    modeling->set_runtime();
}

void Tomography::forward_modeling()
{
    modeling->set_components();

    for (int shot = 0; shot < modeling->total_shots; shot++)
    {
        modeling->shot_id = shot;

        modeling->info_message();

        modeling->initial_setup();

        modeling->forward_solver();

        modeling->build_outputs();

        fill_calculated_data();

        adjoint_state_solver();



    }

    modeling->free_space(); 
}

void Tomography::fill_calculated_data()
{
    for (int node = 0; node < modeling->total_nodes; node++)
    {
        int index = node + modeling->shot_id * modeling->total_nodes;

        dcal[index] = modeling->receiver_output[node]; 
    }
}

void Tomography::adjoint_state_solver()
{
    float cell_volume = modeling->dx * modeling->dy * modeling->dz;

    for (int index = 0; index < n_model; index++) 
    {
        source[index] = 0.0f;    
        adjoint[index] = 1e6f;

        int k = (int) (index / (modeling->nx*modeling->nz));        
        int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;    
        int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz);  

        if ((i == 0) || (i == modeling->nz-1) || 
            (j == 0) || (j == modeling->nx-1) || 
            (k == 0) || (k == modeling->ny-1))  
        {    
            adjoint[index] = 0.0f;        
        }
    }

    for (int node = 0; node < modeling->total_nodes; node++)
    {
        int node_id = node + modeling->shot_id * modeling->total_nodes;

        float tmp = dcal[node_id] + modeling->t0 - dobs[node_id]; 

        int i = (int) floorf(modeling->geometry->nodes.z[node] / modeling->dz);
        int j = (int) floorf(modeling->geometry->nodes.x[node] / modeling->dx);
        int k = (int) floorf(modeling->geometry->nodes.y[node] / modeling->dy);

        int index = i + j*modeling->nz + k*modeling->nz*modeling->nx;

        source[index] += tmp / cell_volume; 
        source[index + 1] += tmp / cell_volume; 
        source[index + modeling->nz] += tmp / cell_volume; 
        source[index + 1 + modeling->nz] += tmp / cell_volume; 
        source[index + modeling->nx * modeling->nz] += tmp / cell_volume; 
        source[index + 1 + modeling->nx * modeling->nz] += tmp / cell_volume; 
        source[index + modeling->nz + modeling->nx * modeling->nz] += tmp / cell_volume; 
        source[index + 1 + modeling->nz + modeling->nx * modeling->nz] += tmp / cell_volume; 
    }





    // for (int sweep = 0; sweep < nSweeps; sweep++)
	// { 
	// 	int start = (sweep == 3 || sweep == 5 || sweep == 6 || sweep == 7) ? totalLevels : meshDim;
	// 	int end = (start == meshDim) ? totalLevels + 1 : meshDim - 1;
	// 	int incr = (start == meshDim) ? true : false;

	// 	int xSweepOff = (sweep == 3 || sweep == 4) ? nxx : 0;
	// 	int ySweepOff = (sweep == 2 || sweep == 5) ? nyy : 0;
	// 	int zSweepOff = (sweep == 1 || sweep == 6) ? nzz : 0;
		
	// 	for (int level = start; level != end; level = (incr) ? level + 1 : level - 1)
	// 	{			
	// 		int xs = max(1, level - (nyy + nzz));	
	// 		int ys = max(1, level - (nxx + nzz));

	// 		int xe = min(nxx, level - (meshDim - 1));
	// 		int ye = min(nyy, level - (meshDim - 1));	
		
	// 		int xr = xe - xs + 1;
	// 		int yr = ye - ys + 1;

	// 		int nThreads = xr * yr;
				
	// 		dim3 bs(16, 16, 1);

	// 		if (nThreads < threadsPerBlock) { bs.x = xr; bs.y = yr; } 

	// 		dim3 gs(iDivUp(xr, bs.x), iDivUp(yr , bs.y), 1);
			
    //         int sgni = sweep + 0*nSweeps;
    //         int sgnj = sweep + 1*nSweeps;
    //         int sgnk = sweep + 2*nSweeps;





	// 	}
	// }









    for (int index = 0; index < n_model; index++) 
    {
        float slowness = 1.0f / modeling->V[index];

        gradient[index] += source[index];
    }
}

void Tomography::compute_residuals()
{
    float r = 0.0f;
    for (int i = 0; i < n_data; i++)
        r += powf(dobs[i] - dcal[i], 2.0f);

    residuo.push_back(sqrtf(r));
}

void Tomography::check_convergence()
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

void Tomography::optimization()
{


}

void Tomography::model_update()
{



}

