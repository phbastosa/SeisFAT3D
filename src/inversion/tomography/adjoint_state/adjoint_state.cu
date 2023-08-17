# include "adjoint_state.cuh"

void Adjoint_State::set_parameters()
{
    set_general_parameters();
    
    set_forward_modeling();

    set_main_components();

    nSweeps = 8;
    meshDim = 3;

	totalLevels = (modeling->nxx - 1) + (modeling->nyy - 1) + (modeling->nzz - 1);

    inversion_method = "[1] - Adjoint State first arrival tomography";

    source = new float[modeling->volsize]();
    adjoint = new float[modeling->volsize]();

    cudaMalloc((void**)&(d_source), modeling->volsize*sizeof(float));
    cudaMalloc((void**)&(d_adjoint), modeling->volsize*sizeof(float));

    std::vector<std::vector<int>> sgnv = {{1,1,1}, {0,1,1}, {1,1,0}, {0,1,0}, {1,0,1}, {0,0,1}, {1,0,0}, {0,0,0}};
    std::vector<std::vector<int>> sgnt = {{1,1,1}, {-1,1,1}, {1,1,-1}, {-1,1,-1}, {1,-1,1}, {-1,-1,1}, {1,-1,-1}, {-1,-1,-1}};

	int * h_sgnv = new int [nSweeps * meshDim]();
	int * h_sgnt = new int [nSweeps * meshDim](); 

	for (int index = 0; index < nSweeps * meshDim; index++)
	{
		int j = index / nSweeps;
		int i = index % nSweeps;				

		h_sgnv[i + j * nSweeps] = sgnv[i][j];
		h_sgnt[i + j * nSweeps] = sgnt[i][j];
	}

	cudaMalloc((void**)&(d_sgnv), nSweeps*meshDim*sizeof(int));
	cudaMalloc((void**)&(d_sgnt), nSweeps*meshDim*sizeof(int));

	cudaMemcpy(d_sgnv, h_sgnv, nSweeps*meshDim*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sgnt, h_sgnt, nSweeps*meshDim*sizeof(int), cudaMemcpyHostToDevice);

    delete[] h_sgnt;
    delete[] h_sgnv;

    std::vector<std::vector<int>>().swap(sgnv);
    std::vector<std::vector<int>>().swap(sgnt);
}

void Adjoint_State::forward_modeling()
{
    initial_setup();

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

void Adjoint_State::initial_setup()
{
    for (int index = 0; index < modeling->nPoints; index++)
    {    
        int k = (int) (index / (modeling->nx*modeling->nz));        
        int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;    
        int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz);          

        int indB = (i+modeling->nbzu) + (j+modeling->nbxl)*modeling->nzz + (k+modeling->nbyl)*modeling->nxx*modeling->nzz;

        modeling->S[indB] = model[index];

        gradient[index] = 0.0f;
    }

    for (int i = 0; i < n_data; i++) dcal[i] = 0.0f; 
}

void Adjoint_State::adjoint_state_solver()
{
    float cell_volume = modeling->dx * modeling->dy * modeling->dz;

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
    }

    for (int node = 0; node < modeling->total_nodes; node++)
    {
        int node_id = node + modeling->shot_id*modeling->total_nodes;

        float tmp = dcal[node_id] + modeling->t0 - dobs[node_id]; 

        int i = (int)(modeling->geometry->nodes.z[node] / modeling->dz);
        int j = (int)(modeling->geometry->nodes.x[node] / modeling->dx);
        int k = (int)(modeling->geometry->nodes.y[node] / modeling->dy);

        int index = i + j*modeling->nzz + k*modeling->nzz*modeling->nxx;

        source[index] += tmp / cell_volume; 
        source[index + 1] += tmp / cell_volume; 
        source[index + modeling->nzz] += tmp / cell_volume; 
        source[index + 1 + modeling->nzz] += tmp / cell_volume; 
        source[index + modeling->nxx*modeling->nzz] += tmp / cell_volume; 
        source[index + 1 + modeling->nxx*modeling->nzz] += tmp / cell_volume; 
        source[index + modeling->nzz + modeling->nxx*modeling->nzz] += tmp / cell_volume; 
        source[index + 1 + modeling->nzz + modeling->nxx*modeling->nzz] += tmp / cell_volume; 
    }

    modeling->expand_boundary(modeling->wavefield_output, modeling->T);

	cudaMemcpy(d_T, modeling->T, modeling->volsize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_source, source, modeling->volsize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_adjoint, adjoint, modeling->volsize*sizeof(float), cudaMemcpyHostToDevice);

    for (int sweep = 0; sweep < nSweeps; sweep++)
	{ 
		int start = (sweep == 3 || sweep == 5 || sweep == 6 || sweep == 7) ? totalLevels : meshDim;
		int end = (start == meshDim) ? totalLevels + 1 : meshDim - 1;
		int incr = (start == meshDim) ? true : false;

		int xSweepOff = (sweep == 3 || sweep == 4) ? modeling->nxx : 0;
		int ySweepOff = (sweep == 2 || sweep == 5) ? modeling->nyy : 0;
		int zSweepOff = (sweep == 1 || sweep == 6) ? modeling->nzz : 0;
		
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

    for (int index = 0; index < modeling->nPoints; index++) 
    {
        int k = (int) (index / (modeling->nx*modeling->nz));        
        int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;    
        int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz);  

        int indp = (i+modeling->nbzu) + (j+modeling->nbxl)*modeling->nzz + (k+modeling->nbyl)*modeling->nxx*modeling->nzz;

        gradient[index] += adjoint[indp]*modeling->S[indp]*modeling->S[indp]*cell_volume / modeling->total_shots;
    }
}

int Adjoint_State::iDivUp(int a, int b) 
{ 
    return ( (a % b) != 0 ) ? (a / b + 1) : (a / b); 
}

void Adjoint_State::optimization()
{







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
		
		if ((z > 0) && (z <= nzz))	
		{
			int i = abs(z - zSweepOffset);
			int j = abs(x - xSweepOffset);
			int k = abs(y - ySweepOffset);

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

                if (d < 1e-6f)
                {
                    adjoint[i + j*nzz + k*nxx*nzz] = 0.0f;    
                }
                else
                {
                    float e = (ap1*adjoint[i + (j-1)*nzz + k*nxx*nzz] - am2*adjoint[i + (j+1)*nzz + k*nxx*nzz]) / dx +
                              (bp1*adjoint[i + j*nzz + (k-1)*nxx*nzz] - bm2*adjoint[i + j*nzz + (k+1)*nxx*nzz]) / dy +
                              (cp1*adjoint[(i-1) + j*nzz + k*nxx*nzz] - cm2*adjoint[(i+1) + j*nzz + k*nxx*nzz]) / dz;

                    float f = (e + source[i + j*nzz + k*nxx*nzz]) / d;
            
                    adjoint[i + j*nzz + k*nxx*nzz] = min(adjoint[i + j*nzz + k*nxx*nzz], f);
                }
            }
        }
    }
}