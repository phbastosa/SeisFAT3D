# include "hfreq_inversion.cuh"

int hfreq_Inversion::iDivUp(int a, int b) 
{ 
    return ( (a % b) != 0 ) ? (a / b + 1) : (a / b); 
}


void hfreq_Inversion::set_forward_modeling()
{
    modeling = new hfreq_Modeling();
}

void hfreq_Inversion::set_inversion_volumes()
{
    n_data = modeling->total_shots * modeling->total_nodes;

    nSweeps = 8;
    meshDim = 3;

    totalLevels = (modeling->nxx - 1) + (modeling->nyy - 1) + (modeling->nzz - 1);

    source = new float[modeling->volsize]();
    adjoint = new float[modeling->volsize]();

    dcal = new float[n_data]();
    dobs = new float[n_data]();

    gradient = new float[modeling->nPoints]();
    slowness = new float[modeling->nPoints]();
    variation = new float[modeling->nPoints]();

    cudaMalloc((void**)&(d_T), modeling->volsize*sizeof(float));
    cudaMalloc((void**)&(d_source), modeling->volsize*sizeof(float));
    cudaMalloc((void**)&(d_adjoint), modeling->volsize*sizeof(float));
}

void hfreq_Inversion::import_obs_data()
{
    int ptr = 0; 
    
    float * data = new float[modeling->total_nodes]();

    for (int shot = 0; shot < modeling->total_shots; shot++)
    {
        import_binary_float(obs_data_folder + obs_data_prefix + std::to_string(shot+1) + ".bin", data, modeling->total_nodes);

        for (int d = ptr; d < ptr + modeling->total_nodes; d++) 
            dobs[d] = data[d - ptr];

        ptr += modeling->total_nodes;        
    }

    delete[] data;        
}

void hfreq_Inversion::get_objective_function()
{
    float square_difference = 0.0f;

    for (int i = 0; i < n_data; i++)
        square_difference += powf(dobs[i] - dcal[i], 2.0f);

    residuo.push_back(sqrtf(square_difference));
}

void hfreq_Inversion::extract_calculated_data()
{
    int skipped = modeling->shot_index * modeling->total_nodes;

    for (int i = 0; i < modeling->total_nodes; i++) 
        dcal[i + skipped] = modeling->receiver_output[i];
}

void hfreq_Inversion::adjoint_propagation()
{
    cell_volume = modeling->dx * modeling->dy * modeling->dz;

    int sidx = (int)(modeling->geometry->shots.x[modeling->shot_index] / modeling->dx) + modeling->nbxl;
    int sidy = (int)(modeling->geometry->shots.y[modeling->shot_index] / modeling->dy) + modeling->nbyl;

    for (int index = 0; index < modeling->volsize; index++) 
    {
        source[index] = 0.0f;    
        adjoint[index] = 1e6f;

        int k = (int) (index / (modeling->nxx*modeling->nzz));        
        int j = (int) (index - k*modeling->nxx*modeling->nzz) / modeling->nzz;    
        int i = (int) (index - j*modeling->nzz - k*modeling->nxx*modeling->nzz);  

        if ((i == 0) || (i == modeling->nzz - 1) || 
            (j == 0) || (j == modeling->nxx - 1) || 
            (k == 0) || (k == modeling->nyy - 1))  
        {    
            adjoint[index] = 0.0f;        
        }
    }

    for (int node = 0; node < modeling->total_nodes; node++)
    {
        int current_node = node + modeling->shot_index*modeling->total_nodes;

        int i = (int)(modeling->geometry->nodes.z[node] / modeling->dz) + modeling->nbzu;
        int j = (int)(modeling->geometry->nodes.x[node] / modeling->dx) + modeling->nbxl;
        int k = (int)(modeling->geometry->nodes.y[node] / modeling->dy) + modeling->nbyl;

        int index = i + j*modeling->nzz + k*modeling->nxx*modeling->nzz;

        if (fabsf(dobs[current_node] - dcal[current_node]) > 1e-3f)
        {
            source[index] += (dobs[current_node] - modeling->T[index]) / cell_volume; 
            
            source[index + 1] += (dobs[current_node] - modeling->T[index + 1]) / cell_volume; 
            source[index - 1] += (dobs[current_node] - modeling->T[index - 1]) / cell_volume; 
            
            source[index + modeling->nzz] += (dobs[current_node] - modeling->T[index + modeling->nzz]) / cell_volume;         
            source[index - modeling->nzz] += (dobs[current_node] - modeling->T[index - modeling->nzz]) / cell_volume;         
            
            source[index + modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index + modeling->nxx*modeling->nzz]) / cell_volume; 
            source[index - modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index - modeling->nxx*modeling->nzz]) / cell_volume; 

            source[index + 1 + modeling->nzz] += (dobs[current_node] - modeling->T[index + 1 + modeling->nzz]) / cell_volume; 
            source[index + 1 - modeling->nzz] += (dobs[current_node] - modeling->T[index + 1 - modeling->nzz]) / cell_volume; 
            source[index - 1 + modeling->nzz] += (dobs[current_node] - modeling->T[index - 1 + modeling->nzz]) / cell_volume; 
            source[index - 1 - modeling->nzz] += (dobs[current_node] - modeling->T[index - 1 - modeling->nzz]) / cell_volume; 
            
            source[index + 1 + modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index + 1 + modeling->nxx*modeling->nzz]) / cell_volume; 
            source[index + 1 - modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index + 1 - modeling->nxx*modeling->nzz]) / cell_volume; 
            source[index - 1 + modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index - 1 + modeling->nxx*modeling->nzz]) / cell_volume; 
            source[index - 1 - modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index - 1 - modeling->nxx*modeling->nzz]) / cell_volume; 
            
            source[index + modeling->nzz + modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index + modeling->nzz + modeling->nxx*modeling->nzz]) / cell_volume; 
            source[index + modeling->nzz - modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index + modeling->nzz - modeling->nxx*modeling->nzz]) / cell_volume; 
            source[index - modeling->nzz + modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index - modeling->nzz + modeling->nxx*modeling->nzz]) / cell_volume; 
            source[index - modeling->nzz - modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index - modeling->nzz - modeling->nxx*modeling->nzz]) / cell_volume; 
            
            source[index + 1 + modeling->nzz + modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index + 1 + modeling->nzz + modeling->nxx*modeling->nzz]) / cell_volume;     
            source[index + 1 + modeling->nzz - modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index + 1 + modeling->nzz - modeling->nxx*modeling->nzz]) / cell_volume;     
            source[index + 1 - modeling->nzz + modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index + 1 - modeling->nzz + modeling->nxx*modeling->nzz]) / cell_volume;     
            source[index + 1 - modeling->nzz - modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index + 1 - modeling->nzz - modeling->nxx*modeling->nzz]) / cell_volume;     
            source[index - 1 + modeling->nzz + modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index - 1 + modeling->nzz + modeling->nxx*modeling->nzz]) / cell_volume;     
            source[index - 1 + modeling->nzz - modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index - 1 + modeling->nzz - modeling->nxx*modeling->nzz]) / cell_volume;     
            source[index - 1 - modeling->nzz + modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index - 1 - modeling->nzz + modeling->nxx*modeling->nzz]) / cell_volume;     
            source[index - 1 - modeling->nzz - modeling->nxx*modeling->nzz] += (dobs[current_node] - modeling->T[index - 1 - modeling->nzz - modeling->nxx*modeling->nzz]) / cell_volume;     
        }
    }

    cudaMemcpy(d_T, modeling->T, modeling->volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_source, source, modeling->volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_adjoint, adjoint, modeling->volsize*sizeof(float), cudaMemcpyHostToDevice);

    // for (int sweepCount = 0; sweepCount < meshDim; sweepCount++)
    // {
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
    // }
    
    cudaMemcpy(adjoint, d_adjoint, modeling->volsize*sizeof(float), cudaMemcpyDeviceToHost);

    for (int index = 0; index < modeling->nPoints; index++) 
    {
        int k = (int) (index / (modeling->nx*modeling->nz));        
        int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;    
        int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz);  

        int indp = (i + modeling->nbzu) + (j + modeling->nbxl)*modeling->nzz + (k + modeling->nbyl)*modeling->nxx*modeling->nzz;

        gradient[index] += adjoint[indp]*modeling->T[indp]*modeling->S[indp]*modeling->S[indp]*cell_volume / modeling->total_shots / modeling->total_nodes;
    }    
}

void hfreq_Inversion::update_specifications()
{
    for (int index = 0; index < modeling->nPoints; index++)
    {
        slowness[index] += variation[index];

        modeling->model[index] = 1.0f / slowness[index];
    }

    modeling->expand_boundary(slowness, modeling->S);
}

__global__ void adjoint_state_kernel(float * adjoint, float * source, float * T, int level, int xOffset, int yOffset, int xSweepOffset, int ySweepOffset, 
                                     int zSweepOffset, int nxx, int nyy, int nzz, float dx, float dy, float dz)
{
    int x = (blockIdx.x * blockDim.x + threadIdx.x) + xOffset;
    int y = (blockIdx.y * blockDim.y + threadIdx.y) + yOffset;

    if ((x < nxx) && (y < nyy)) 
    {
        int z = level - (x + y);
		
	if ((z >= 0) && (z < nzz))	
	{
	    int i = (int)abs(z - zSweepOffset);
	    int j = (int)abs(x - xSweepOffset);
	    int k = (int)abs(y - ySweepOffset);

            if ((i > 0) && (i < nzz - 1) && (j > 0) && (j < nxx - 1) && (k > 0) && (k < nyy - 1))
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

                if (fabsf(d) < 1e-6f)
                {
                    adjoint[i + j*nzz + k*nxx*nzz] = 0.0f;    
                }
                else
                {
                    float e = (ap1*adjoint[i + (j-1)*nzz + k*nxx*nzz] - am2*adjoint[i + (j+1)*nzz + k*nxx*nzz]) / dx +
                              (bp1*adjoint[i + j*nzz + (k-1)*nxx*nzz] - bm2*adjoint[i + j*nzz + (k+1)*nxx*nzz]) / dy +
                              (cp1*adjoint[(i-1) + j*nzz + k*nxx*nzz] - cm2*adjoint[(i+1) + j*nzz + k*nxx*nzz]) / dz;

                    float f = (e + source[i + j*nzz + k*nxx*nzz]) / d;

                    if (adjoint[i + j*nzz + k*nxx*nzz] > f) 
                        adjoint[i + j*nzz + k*nxx*nzz] = f;
                }
            }
        }
    }
}
