# include "adjoint_state.cuh"

int Adjoint_State::iDivUp(int a, int b) 
{ 
    return ( (a % b) != 0 ) ? (a / b + 1) : (a / b); 
}

void Adjoint_State::set_parameters()
{
    set_general_parameters();
    
    set_forward_modeling();

    set_main_components();

    nSweeps = 8;
    meshDim = 3;

    max_slowness_variation = 1e-4f;

	totalLevels = (modeling->nxx - 1) + (modeling->nyy - 1) + (modeling->nzz - 1);

    inversion_method = "[1] - Adjoint State first arrival tomography";

    source = new float[modeling->volsize]();
    adjoint = new float[modeling->volsize]();

    cudaMalloc((void**)&(d_T), modeling->volsize*sizeof(float));
    cudaMalloc((void**)&(d_source), modeling->volsize*sizeof(float));
    cudaMalloc((void**)&(d_adjoint), modeling->volsize*sizeof(float));
}

void Adjoint_State::forward_modeling()
{
    adjoint_initial_setup();

    for (int shot = 0; shot < modeling->total_shots; shot++)
    {
        modeling->shot_id = shot;
    
        modeling->info_message();

        tomography_message();

        modeling->initial_setup();
        modeling->forward_solver();
        modeling->build_outputs();

        extract_calculated_data();
        
        adjoint_state_solver();
    }

    export_gradient();
}

void Adjoint_State::adjoint_initial_setup()
{
    # pragma omp parallel for
    for (int index = 0; index < modeling->nPoints; index++) 
    {
        int k = (int) (index / (modeling->nx*modeling->nz));        
        int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;    
        int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz);          

        int indB = (i+modeling->nbzu) + (j+modeling->nbxl)*modeling->nzz + (k+modeling->nbyl)*modeling->nxx*modeling->nzz;

        modeling->S[indB] = model[index];

        gradient[index] = 0.0f;
    }

    # pragma omp parallel for
    for (int i = 0; i < n_data; i++) 
        dcal[i] = 0.0f; 
}

void Adjoint_State::adjoint_state_solver()
{
    Tmax = 0.0f;

    cell_volume = modeling->dx * modeling->dy * modeling->dz;

    modeling->expand_boundary(modeling->wavefield_output, modeling->T);

    # pragma omp parallel for
    for (int index = 0; index < modeling->volsize; index++) 
    {
        source[index] = 0.0f;    
        adjoint[index] = 1e6f;

        int k = (int) (index / (modeling->nxx*modeling->nzz));        
        int j = (int) (index - k*modeling->nxx*modeling->nzz) / modeling->nzz;    
        int i = (int) (index - j*modeling->nzz - k*modeling->nxx*modeling->nzz);  

        if ((i == 0) || (i == modeling->nzz-1) || 
            (j == 0) || (j == modeling->nxx-1) || 
            (k == 0) || (k == modeling->nyy-1))  
        {    
            adjoint[index] = 0.0f;        
        }

        if (Tmax < modeling->T[index]) Tmax = modeling->T[index];    
    }

    for (int node = 0; node < modeling->total_nodes; node++)
    {
        int current_node = node + modeling->shot_id*modeling->total_nodes;

        int i = (int)(modeling->geometry->nodes.z[node] / modeling->dz) + modeling->nbzu;
        int j = (int)(modeling->geometry->nodes.x[node] / modeling->dx) + modeling->nbxl;
        int k = (int)(modeling->geometry->nodes.y[node] / modeling->dy) + modeling->nbyl;

        int index = i + j*modeling->nzz + k*modeling->nxx*modeling->nzz;

        source[index] += (dobs[current_node] - modeling->T[index]) / cell_volume; 
        source[index + 1] += (dobs[current_node] - modeling->T[index + 1]) / cell_volume; 
        source[index + modeling->nzz] += (dobs[current_node] - modeling->T[index + modeling->nzz]) / cell_volume;         
        source[index + 1 + modeling->nzz] += (dobs[current_node] - modeling->T[index + 1 + modeling->nzz]) / cell_volume; 
        source[index + modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index + modeling->nxx*modeling->nzz]) / cell_volume; 
        source[index + 1 + modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index + 1 + modeling->nxx*modeling->nzz]) / cell_volume; 
        source[index + modeling->nzz + modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index + modeling->nzz + modeling->nxx*modeling->nzz]) / cell_volume; 
        source[index + 1 + modeling->nzz + modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index + 1 + modeling->nzz + modeling->nxx*modeling->nzz]) / cell_volume;     
    }

	cudaMemcpy(d_T, modeling->T, modeling->volsize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_source, source, modeling->volsize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_adjoint, adjoint, modeling->volsize*sizeof(float), cudaMemcpyHostToDevice);

    for (int sweep = 0; sweep < nSweeps; sweep++)
	{ 
		int start = (sweep == 3 || sweep == 5 || sweep == 6 || sweep == 7) ? totalLevels : meshDim;
		int end = (start == meshDim) ? totalLevels + 1 : meshDim - 1;
		int incr = (start == meshDim) ? true : false;

		int xSweepOff = (sweep == 3 || sweep == 4) ? modeling->nxx + 1 : 0;
		int ySweepOff = (sweep == 2 || sweep == 5) ? modeling->nyy + 1 : 0;
		int zSweepOff = (sweep == 1 || sweep == 6) ? modeling->nzz + 1 : 0;
		
		for (int level = start; level != end; level = (incr) ? level + 1 : level - 1)
		{			
			int xs = max(1, level - (modeling->nyy + modeling->nzz));	
			int ys = max(1, level - (modeling->nxx + modeling->nzz));

			int xe = min(modeling->nxx, level - (meshDim - 1));
			int ye = min(modeling->nyy, level - (meshDim - 1));	
		
			int xr = xe - xs + 1;
			int yr = ye - ys + 1;

			int nThreads = xr * yr;
				
			dim3 bs(16, 16, 1);

			if (nThreads < modeling->threadsPerBlock) { bs.x = xr; bs.y = yr; } 

			dim3 gs(iDivUp(xr, bs.x), iDivUp(yr , bs.y), 1);

            adjoint_state_kernel<<<gs,bs>>>(d_adjoint, d_source, d_T, level, xs, ys, xSweepOff, ySweepOff, zSweepOff, 
                                            modeling->nxx, modeling->nyy, modeling->nzz, modeling->dx, modeling->dy, modeling->dz);

            cudaDeviceSynchronize();
		}
	}

    cudaMemcpy(adjoint, d_adjoint, modeling->volsize*sizeof(float), cudaMemcpyDeviceToHost);

    adjoint_conditioning();

    # pragma omp parallel for
    for (int index = 0; index < modeling->nPoints; index++) 
    {
        int k = (int) (index / (modeling->nx*modeling->nz));        
        int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;    
        int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz);  

        int indp = (i + modeling->nbzu) + (j + modeling->nbxl)*modeling->nzz + (k + modeling->nbyl)*modeling->nxx*modeling->nzz;

        gradient[index] += adjoint[indp]*modeling->S[indp]*modeling->S[indp]*cell_volume;
    }
}

void Adjoint_State::adjoint_conditioning() 
{
    int sidx = (int)(modeling->geometry->shots.x[modeling->shot_id] / modeling->dx) + modeling->nbxl;
    int sidy = (int)(modeling->geometry->shots.y[modeling->shot_id] / modeling->dy) + modeling->nbyl;

    # pragma omp parallel for
    for (int i = 0; i < modeling->nzz; i++)
    {
        for (int j = 0; j < modeling->nxx; j++)
        {
            adjoint[i + j*modeling->nzz + sidy*modeling->nxx*modeling->nzz] = 0.5f*(adjoint[i + j*modeling->nzz + (sidy+1)*modeling->nxx*modeling->nzz] + adjoint[i + j*modeling->nzz + (sidy-1)*modeling->nxx*modeling->nzz]);
        }

        for (int k = 0; k < modeling->nyy; k++)    
        {
            adjoint[i + sidx*modeling->nzz + k*modeling->nxx*modeling->nzz] = 0.5f*(adjoint[i + (sidx+1)*modeling->nzz + k*modeling->nxx*modeling->nzz] + adjoint[i + (sidx-1)*modeling->nzz + k*modeling->nxx*modeling->nzz]);
        }
    }
    
    # pragma omp parallel for
    for (int index = 0; index < modeling->volsize; index++)
    {
        adjoint[index] *= modeling->T[index] / Tmax;
    }
}

void Adjoint_State::optimization() 
{ 
    gradient_normalization();
  
    backtracking_linesearch();

    nonlinear_conjugate_gradient();
}

void Adjoint_State::gradient_normalization()
{
    norm = 0.0f;

    for (int index = 0; index < modeling->nPoints; index++)
    {
        if (norm < abs(gradient[index]))
            norm = abs(gradient[index]);        
    }    

    norm /= max_slowness_variation;
}

void Adjoint_State::backtracking_linesearch()
{
    std::cout<<"\nRunning backtracking linesearch\n"<<std::endl;

    float k1 = 1e-4f;
    float k2 = 0.85f;

    float dot_product = 0.0f;

    # pragma omp parallel for reduction(+:dot_product)
    for (int index = 0; index < modeling->nPoints; index++)
        dot_product += (gradient[index]) * (gradient[index]);       

    alpha = 1.0f;

    float f0 = residuo.back();

    float f1 = get_objective_function(alpha, gradient);

    std::cout<<f1<<std::endl;
    std::cout<<f0 - k1*alpha*dot_product<<std::endl;

    //while (f1 > f0 + k1*alpha*dot_product)
    //alpha *= k2;
}

float Adjoint_State::get_objective_function(float step, float * grad)
{
    # pragma omp parallel for
    for (int index = 0; index < modeling->nPoints; index++)
    {
        int k = (int) (index / (modeling->nx*modeling->nz));        
        int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;    
        int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz);          

        int indB = (i + modeling->nbzu) + (j + modeling->nbxl)*modeling->nzz + (k + modeling->nbyl)*modeling->nxx*modeling->nzz;
        
        modeling->S[indB] = model[index] + step*grad[index]/norm;
    }

    for (int shot = 0; shot < modeling->total_shots; shot++)
    {
        modeling->shot_id = shot;

        modeling->initial_setup();
        modeling->forward_solver();
        modeling->build_outputs();

        extract_calculated_data();
    }

    float function = 0.0f;

    for (int i = 0; i < n_data; i++)
        function += powf(dobs[i] - dcal[i], 2.0f);

    function = sqrtf(function);

    return function;
}

void Adjoint_State::nonlinear_conjugate_gradient()
{
    // if (iteration == 1)

    // else

    # pragma omp parallel for    
    for (int index = 0; index < modeling->nPoints; index++)
        dm[index] = alpha*gradient[index]/norm;





}

__global__ void adjoint_state_kernel(float * adjoint, float * source, float * T, int level, int xOffset, int yOffset, 
                                     int xSweepOffset, int ySweepOffset, int zSweepOffset, int nxx, int nyy, int nzz, 
                                     float dx, float dy, float dz)
{
	int x = (blockIdx.x * blockDim.x + threadIdx.x) + xOffset;
	int y = (blockIdx.y * blockDim.y + threadIdx.y) + yOffset;

	if ((x <= nxx) && (y <= nyy)) 
	{
		int z = level - (x + y);
		
		if ((z >= 0) && (z <= nzz))	
		{
			int i = (int)abs(z - zSweepOffset);
			int j = (int)abs(x - xSweepOffset);
			int k = (int)abs(y - ySweepOffset);

			if ((i > 0) && (i < nzz-1) && (j > 0) && (j < nxx-1) && (k > 0) && (k < nyy-1))
			{		
                float a1 = -1.0f*(T[i + j*nzz + k*nxx*nzz] - T[i + (j-1)*nzz + k*nxx*nzz]) / dx;
                float ap1 = (a1 + abs(a1)) / 2.0f;
                float am1 = (a1 - abs(a1)) / 2.0f;

                float a2 = -1.0f*(T[i + (j+1)*nzz + k*nxx*nzz] - T[i + j*nzz + k*nxx*nzz]) / dx;
                float ap2 = (a2 + abs(a2)) / 2.0f;
                float am2 = (a2 - abs(a2)) / 2.0f;

                float b1 = -1.0f*(T[i + j*nzz + k*nxx*nzz] - T[i + j*nzz + (k-1)*nxx*nzz]) / dy;
                float bp1 = (b1 + abs(b1)) / 2.0f;
                float bm1 = (b1 - abs(b1)) / 2.0f;

                float b2 = -1.0f*(T[i + j*nzz + (k+1)*nxx*nzz] - T[i + j*nzz + k*nxx*nzz]) / dy;
                float bp2 = (b2 + abs(b2)) / 2.0f;
                float bm2 = (b2 - abs(b2)) / 2.0f;

                float c1 = -1.0f*(T[i + j*nzz + k*nxx*nzz] - T[(i-1) + j*nzz + k*nxx*nzz]) / dz;
                float cp1 = (c1 + abs(c1)) / 2.0f;
                float cm1 = (c1 - abs(c1)) / 2.0f;

                float c2 = -1.0f*(T[(i+1) + j*nzz + k*nxx*nzz] - T[i + j*nzz + k*nxx*nzz]) / dz;
                float cp2 = (c2 + abs(c2)) / 2.0f;
                float cm2 = (c2 - abs(c2)) / 2.0f;

                float d = (ap2 - am1)/dx + (bp2 - bm1)/dy + (cp2 - cm1)/dz;

                if (abs(d) < 1e-6f)
                {
                    adjoint[i + j*nzz + k*nxx*nzz] = 0.0f;    
                }
                else
                {
                    float e = (ap1*adjoint[i + (j-1)*nzz + k*nxx*nzz] - am2*adjoint[i + (j+1)*nzz + k*nxx*nzz]) / dx +
                              (bp1*adjoint[i + j*nzz + (k-1)*nxx*nzz] - bm2*adjoint[i + j*nzz + (k+1)*nxx*nzz]) / dy +
                              (cp1*adjoint[(i-1) + j*nzz + k*nxx*nzz] - cm2*adjoint[(i+1) + j*nzz + k*nxx*nzz]) / dz;

                    float f = (e + source[i + j*nzz + k*nxx*nzz]) / d;
                    float g = adjoint[i + j*nzz + k*nxx*nzz];

                    if (g > f) adjoint[i + j*nzz + k*nxx*nzz] = f;
                }
            }
        }
    }
}