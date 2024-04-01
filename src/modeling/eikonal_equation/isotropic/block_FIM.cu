# include "block_FIM.cuh"

void Block_FIM::set_specifics()
{
    int pdx = (BLOCK_LENGTH - nx % BLOCK_LENGTH) % BLOCK_LENGTH;
    int pdy = (BLOCK_LENGTH - ny % BLOCK_LENGTH) % BLOCK_LENGTH;
    int pdz = (BLOCK_LENGTH - nz % BLOCK_LENGTH) % BLOCK_LENGTH;    
    
    if (pdx == 0)      { nbxl = 2; nbxr = 2; }
    else if (pdx == 1) { nbxl = 3; nbxr = 2; }
    else if (pdx == 2) { nbxl = 3; nbxr = 3; }
    else if (pdx == 3) { nbxl = 4; nbxr = 3; }

    if (pdy == 0)      { nbyl = 2; nbyr = 2; }
    else if (pdy == 1) { nbyl = 3; nbyr = 2; }
    else if (pdy == 2) { nbyl = 3; nbyr = 3; }
    else if (pdy == 3) { nbyl = 4; nbyr = 3; }

    if (pdz == 0)      { nbzu = 2; nbzd = 2; }
    else if (pdz == 1) { nbzu = 3; nbzd = 2; }
    else if (pdz == 2) { nbzu = 3; nbzd = 3; }
    else if (pdz == 3) { nbzu = 4; nbzd = 3; }
}

void Block_FIM::set_volumes()
{
    type_name = std::string("fim");
    type_message = std::string("[1] - Block FIM (Jeong & Whitaker, 2008)");

    blksize = BLOCK_LENGTH * BLOCK_LENGTH * BLOCK_LENGTH;

    nbx = nxx / BLOCK_LENGTH;
    nby = nyy / BLOCK_LENGTH;
    nbz = nzz / BLOCK_LENGTH;
	
    nblk = nbx * nby * nbz;

    nit = 10;

    h_mask = new bool[volsize]();
    h_slow = new float[volsize]();
    h_time = new float[volsize]();
	
    h_list = new uint[nblk]();      
    h_listed = new bool[nblk]();    
    h_listVol = new bool[nblk]();   

    if ((dx != dy) || (dx != dz))
        throw std::invalid_argument("\033[31mError: Block_FIM algorithm must to be the model spacing fixed (dx = dy = dz).\033[0;0m");

    T = new float[volsize](); 

    cudaMalloc((void**)&(d_con), volsize*sizeof(bool));  
    cudaMalloc((void**)&(d_mask), volsize*sizeof(bool));
    cudaMalloc((void**)&(d_slow), volsize*sizeof(float));
    cudaMalloc((void**)&(d_time), volsize*sizeof(float));
    cudaMalloc((void**)&(t_time), volsize*sizeof(float)); 
    cudaMalloc((void**)&(d_list), nblk*sizeof(uint));
    cudaMalloc((void**)&(d_listVol), nblk*sizeof(bool));
}

void Block_FIM::initialization()
{
    uint idx = 0;
    uint blk_idx = 0;
    uint list_idx = 0;
    
    nActiveBlock = 0;

    float t0 = S[source_index] * sqrtf(powf((float)(sidx*dx) - geometry->shots.x[shot_index], 2.0f) +
                                       powf((float)(sidy*dy) - geometry->shots.y[shot_index], 2.0f) +
                                       powf((float)(sidz*dz) - geometry->shots.z[shot_index], 2.0f));

    for(int zStr = 0; zStr < nzz; zStr += BLOCK_LENGTH) 
    {
        for(int yStr = 0; yStr < nyy; yStr += BLOCK_LENGTH) 
	{
      	    for(int xStr = 0; xStr < nxx; xStr += BLOCK_LENGTH) 
	    {
        	bool isSeedBlock = false;

        	for(int z = zStr; z < zStr + BLOCK_LENGTH; z++) 
		{
		    for(int y = yStr; y < yStr + BLOCK_LENGTH; y++) 
		    {
            	        for(int x = xStr; x < xStr + BLOCK_LENGTH; x++) 
			{
                            h_time[idx] = INF;
                            h_mask[idx] = true;
                            h_slow[idx] = S[z + x*nzz + y*nxx*nzz];

                            if (x == sidx && y == sidy && z == sidz) 
                            {
                                h_time[idx] = t0;
                                isSeedBlock = true;
                            }

			    ++idx;
            	        }
          	    }
	        }
        
                if(isSeedBlock) 
	        {          			
          	    h_list[list_idx] = blk_idx;
          	    h_listed[blk_idx] = true;
		    h_listVol[blk_idx] = true;
          			
                    ++list_idx;
                    ++nActiveBlock;
                } 
	        else 
	        {
          	    h_listed[blk_idx] = false;
          	    h_listVol[blk_idx] = false;
                }
        		
	        ++blk_idx;
      	    }
        }	
    }
    
    cudaMemset(d_con, 1, volsize*sizeof(bool));
    cudaMemcpy(d_slow, h_slow, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_time, h_time, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(t_time, h_time, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_mask, h_mask, volsize*sizeof(bool), cudaMemcpyHostToDevice);
    cudaMemcpy(d_listVol, h_listVol, nblk*sizeof(bool), cudaMemcpyHostToDevice);
}

void Block_FIM::set_forward_solver()
{
    float dh = dx;

    int nTotalIter = 0;
    uint nTotalBlockProcessed = 0;

    dim3 dimEntireGrid(nblk);
    dim3 dimGrid(nActiveBlock);
    dim3 dimBlock(BLOCK_LENGTH, BLOCK_LENGTH, BLOCK_LENGTH);

    while(nActiveBlock > 0)
    {
        assert(nActiveBlock < 4294967295);

        nTotalBlockProcessed += nActiveBlock;

        nTotalIter++;

        // 1. Run solver on current active tiles.

        dimGrid.y = (unsigned int)floorf((float)(nActiveBlock-1)/65535)+1;
        dimGrid.x = (unsigned int)ceilf((float)nActiveBlock/(float)dimGrid.y);

        cudaMemcpy(d_list, h_list, nActiveBlock*sizeof(uint), cudaMemcpyHostToDevice);

        run_solver<<<dimGrid,dimBlock>>>(d_slow, d_mask, d_time, t_time, d_con, d_list, nxx, nyy, nzz, dh, nit, nActiveBlock);        
        cudaDeviceSynchronize();

        // 2. Reduction.

        run_reduction<<<dimGrid,dim3(BLOCK_LENGTH,BLOCK_LENGTH,BLOCK_LENGTH/2)>>>(d_con, d_listVol, d_list, nActiveBlock);
        cudaDeviceSynchronize();

        // 3. Check neighbor tiles of converged tile.
        // Add any active block of neighbor of converged block is inserted to the list.        
        
        cudaMemcpy(h_listVol, d_listVol, nblk*sizeof(bool), cudaMemcpyDeviceToHost);

        uint nBlkX = nxx/BLOCK_LENGTH;
        uint nBlkY = nyy/BLOCK_LENGTH;
        uint nOldActiveBlock = nActiveBlock;

        for(uint i = 0; i < nOldActiveBlock; i++)
        {
            // check 6-neighbor of current active tile
            uint currBlkIdx = h_list[i];

            if(!h_listVol[currBlkIdx]) // not active : converged
            {
                uint nb[6];

                nb[0] = (currBlkIdx < nBlkX*nBlkY) ? currBlkIdx : (currBlkIdx - nBlkX*nBlkY);  //tp
                nb[1] = ((currBlkIdx + nBlkX*nBlkY) >= nblk) ? currBlkIdx : (currBlkIdx + nBlkX*nBlkY); //bt
                nb[2] = (currBlkIdx < nBlkX) ? currBlkIdx : (currBlkIdx - nBlkX); //up
                nb[3] = ((currBlkIdx + nBlkX) >= nblk) ? currBlkIdx : (currBlkIdx + nBlkX); //dn
                nb[4] = (currBlkIdx%nBlkX == 0) ? currBlkIdx : currBlkIdx-1; //lf
                nb[5] = ((currBlkIdx+1)%nBlkX == 0) ? currBlkIdx : currBlkIdx+1; //rt

                for(int nbIdx = 0; nbIdx < 6; nbIdx++)
                {
                    uint currIdx = nb[nbIdx];

                    if(!h_listed[currIdx])
                    {
                        h_listed[currIdx] = true;
                        h_list[nActiveBlock++] = currIdx;
                    }
                }
            }
        }

        // 4. Run solver only once for neighbor blocks of converged block.
        // Current active list contains active blocks and neighbor blocks of any converged blocks.

        // Update grid dimension because nActiveBlock is changed
        dimGrid.y = (unsigned int)floor(((float)nActiveBlock-1)/65535)+1;
        dimGrid.x = (unsigned int)ceil((float)nActiveBlock/(float)dimGrid.y);

        cudaMemcpy(d_list, h_list, nActiveBlock*sizeof(uint), cudaMemcpyHostToDevice);
        
        run_check_neighbor<<<dimGrid,dimBlock>>>(d_slow, d_mask, t_time, d_time, d_con, d_list, nxx, nyy, nzz, dh, nOldActiveBlock, nActiveBlock);
        cudaDeviceSynchronize();

        // 5. Reduction.

        run_reduction<<<dimGrid,dim3(BLOCK_LENGTH,BLOCK_LENGTH,BLOCK_LENGTH/2)>>>(d_con, d_listVol, d_list, nActiveBlock);
        cudaDeviceSynchronize();

        // 6. Update active list.
        // Read back active volume from the device and add active block to active list on the host memory.

        nActiveBlock = 0;
        
        cudaMemcpy(h_listVol, d_listVol, nblk*sizeof(bool), cudaMemcpyDeviceToHost);

        for(uint i = 0; i < nblk; i++)
        {
            if(h_listVol[i]) // true : active block (not converged)
            {
                h_listed[i] = true;
                h_list[nActiveBlock++] = i;
            }
            else
            { 
                h_listed[i] = false;
            }
        }

        cudaDeviceSynchronize();
    }

    cudaMemcpy(h_time, d_time, volsize*sizeof(float), cudaMemcpyDeviceToHost);

    uint idx = 0;

    for(int zStr = 0; zStr < nzz; zStr += BLOCK_LENGTH) 
    {
        for(int yStr = 0; yStr < nyy; yStr += BLOCK_LENGTH) 
	{
	    for(int xStr = 0; xStr < nxx; xStr += BLOCK_LENGTH) 
	    {
	        for(int z = zStr; z < zStr + BLOCK_LENGTH; z++) 
		{
		    for(int y = yStr; y < yStr + BLOCK_LENGTH; y++) 
		    {
		        for(int x = xStr; x < xStr + BLOCK_LENGTH; x++) 
			{
                            T[z + x*nzz + y*nxx*nzz] = h_time[++idx];
			}
		    }
		}
	    }
	}
    }

    get_wavefield_output();	
    get_receiver_output();
}

void Block_FIM::free_space()
{
	delete[] h_mask;
	delete[] h_slow;
	delete[] h_time;
	
	delete[] h_list;
	delete[] h_listed;
	delete[] h_listVol;

	cudaFree(d_con);
	cudaFree(d_list);
	cudaFree(d_listVol);

	cudaFree(d_mask);
	cudaFree(d_slow);
	cudaFree(d_time);
	cudaFree(t_time);
}

__global__ void run_solver(float* spd, bool* mask, const float *sol_in, float *sol_out, bool *con, uint* list, int xdim, int ydim, int zdim, float dh, int nIter, uint nActiveBlock)
{
    uint list_idx = blockIdx.y*gridDim.x + blockIdx.x;

    if(list_idx < nActiveBlock)
    {
        // retrieve actual block index from the active list
        uint block_idx = list[list_idx];

	uint blocksize = BLOCK_LENGTH*BLOCK_LENGTH*BLOCK_LENGTH;
	uint base_addr = block_idx*blocksize;

	uint xgridlength = xdim/BLOCK_LENGTH;
	uint ygridlength = ydim/BLOCK_LENGTH;
	uint zgridlength = zdim/BLOCK_LENGTH;

	// compute block index
	uint bx = block_idx % xgridlength;
	uint tmpIdx = (block_idx - bx) / xgridlength;
	uint by = tmpIdx % ygridlength;
	uint bz = (tmpIdx - by) / ygridlength;

	uint tx = threadIdx.x;
	uint ty = threadIdx.y;
	uint tz = threadIdx.z;

	uint tIdx = tz*BLOCK_LENGTH*BLOCK_LENGTH + ty*BLOCK_LENGTH + tx;

	__shared__ float _sol[BLOCK_LENGTH+2][BLOCK_LENGTH+2][BLOCK_LENGTH+2];

	// copy global to shared memory
	dim3 idx(tx+1,ty+1,tz+1);

	SOL(idx.x,idx.y,idx.z) = sol_in[base_addr + tIdx];
		
        float F = spd[base_addr + tIdx];

	bool isValid = mask[base_addr + tIdx];

	uint new_base_addr, new_tIdx;

	// 1-neighborhood values
	if(tx == 0) 
	{
	    if(bx == 0) // end of the grid
	    {	
	        new_tIdx = tIdx;
		new_base_addr = base_addr;
	    }
	    else
	    {
	        new_tIdx = tIdx + BLOCK_LENGTH-1;
		new_base_addr = (block_idx - 1)*blocksize;	
	    }

	    SOL(tx,idx.y,idx.z) = sol_in[new_base_addr + new_tIdx];	
	}

	if(tx == BLOCK_LENGTH-1)
	{
	    if(bx == xgridlength-1) // end of the grid
	    {
	        new_tIdx = tIdx;
		new_base_addr = base_addr;
	    }
	    else
	    {
	        new_tIdx = tIdx - (BLOCK_LENGTH-1);
		new_base_addr = (block_idx + 1)*blocksize;	
	    }
	    
            SOL(tx+2,idx.y,idx.z) = sol_in[new_base_addr + new_tIdx];	
	}

	if(ty == 0)
	{
	    if(by == 0)
	    {
	        new_tIdx = tIdx;
		new_base_addr = base_addr;
	    }
	    else
	    {
	        new_tIdx = tIdx + (BLOCK_LENGTH-1)*BLOCK_LENGTH;
		new_base_addr = (block_idx - xgridlength)*blocksize;
	    }

	    SOL(idx.x,ty,idx.z) = sol_in[new_base_addr + new_tIdx];
	}

	if(ty == BLOCK_LENGTH-1)
	{
	    if(by == ygridlength-1) 
	    {
	        new_tIdx = tIdx;
		new_base_addr = base_addr;
	    }
	    else
	    {
	        new_tIdx = tIdx - (BLOCK_LENGTH-1)*BLOCK_LENGTH;
		new_base_addr = (block_idx + xgridlength)*blocksize;
	    }

	    SOL(idx.x,ty+2,idx.z) = sol_in[new_base_addr + new_tIdx];
	}

	if(tz == 0)
	{
	    if(bz == 0)
	    {
	        new_tIdx = tIdx;
		new_base_addr = base_addr;
	    }
	    else
	    {
	        new_tIdx = tIdx + (BLOCK_LENGTH-1)*BLOCK_LENGTH*BLOCK_LENGTH;
		new_base_addr = (block_idx - xgridlength*ygridlength)*blocksize;
	    }

	    SOL(idx.x,idx.y,tz) = sol_in[new_base_addr + new_tIdx];
	}

	if(tz == BLOCK_LENGTH-1)
	{
	    if(bz == zgridlength-1) 
	    {
	        new_tIdx = tIdx;
		new_base_addr = base_addr;
	    }
	    else
	    {
	        new_tIdx = tIdx - (BLOCK_LENGTH-1)*BLOCK_LENGTH*BLOCK_LENGTH;
		new_base_addr = (block_idx + xgridlength*ygridlength)*blocksize;
	    }

	    SOL(idx.x,idx.y,tz+2) = sol_in[new_base_addr + new_tIdx];
	}

	__syncthreads();

	float a,b,c,oldT,newT;

	for(int iter=0; iter<nIter; iter++)	
	{
	    // compute new value
	    oldT = newT = SOL(idx.x,idx.y,idx.z);

	    if(isValid)
	    {
		a = min(SOL(tx,idx.y,idx.z),SOL(tx+2,idx.y,idx.z));
		b = min(SOL(idx.x,ty,idx.z),SOL(idx.x,ty+2,idx.z));
		c = min(SOL(idx.x,idx.y,tz),SOL(idx.x,idx.y,tz+2));

		float tmp = get_time_eikonal(a, b, c, dh, F);

		newT = min(tmp,oldT);
	    }
	    __syncthreads();	

	    if(isValid) SOL(idx.x,idx.y,idx.z) = newT;
        }

	float residue = oldT - newT;

	// write back to global memory
	con[base_addr + tIdx] = (residue < EPS) ? true : false;
	sol_out[base_addr + tIdx] = newT;		
    }
}

__global__ void run_reduction(bool *con, bool *listVol, uint *list, uint nActiveBlock)
{
    uint list_idx = blockIdx.y*gridDim.x + blockIdx.x;

    if(list_idx < nActiveBlock)
    {
        uint block_idx = list[list_idx];

	__shared__ bool conv[BLOCK_LENGTH*BLOCK_LENGTH*BLOCK_LENGTH];

	uint blocksize = BLOCK_LENGTH*BLOCK_LENGTH*BLOCK_LENGTH/2;
	uint base_addr = block_idx*blocksize*2;
		
        uint tx = threadIdx.x;
	uint ty = threadIdx.y;
	uint tz = threadIdx.z;

	uint tIdx = tz*BLOCK_LENGTH*BLOCK_LENGTH + ty*BLOCK_LENGTH + tx;

	conv[tIdx] = con[base_addr + tIdx];
	conv[tIdx + blocksize] = con[base_addr + tIdx + blocksize];

	__syncthreads();

	for(uint i=blocksize; i>0; i/=2)
	{
	    if(tIdx < i)
	    {
		bool b1, b2;
		b1 = conv[tIdx];
		b2 = conv[tIdx+i];
		conv[tIdx] = (b1 && b2) ? true : false ;
	    }
	    __syncthreads();
	}

        // active list is negation of tile convergence (active = not converged)
        if(tIdx == 0) listVol[block_idx] = !conv[0]; 
    }
}

__global__ void run_check_neighbor(float* spd, bool* mask, const float *sol_in, float *sol_out, bool *con, uint* list, int xdim, int ydim, int zdim, float dh, uint nActiveBlock, uint nTotalBlock)
{
    uint list_idx = blockIdx.y*gridDim.x + blockIdx.x;

    if(list_idx < nTotalBlock)
    {
        __shared__ float _sol[BLOCK_LENGTH+2][BLOCK_LENGTH+2][BLOCK_LENGTH+2];

	uint block_idx = list[list_idx];
	uint blocksize = BLOCK_LENGTH*BLOCK_LENGTH*BLOCK_LENGTH;
	uint base_addr = block_idx*blocksize;

	uint tx = threadIdx.x;
	uint ty = threadIdx.y;
	uint tz = threadIdx.z;
	uint tIdx = tz*BLOCK_LENGTH*BLOCK_LENGTH + ty*BLOCK_LENGTH + tx;

	if(list_idx < nActiveBlock) // copy value
	{
	    sol_out[base_addr + tIdx] = sol_in[base_addr + tIdx];
	} 
	else
	{
	    uint xgridlength = xdim/BLOCK_LENGTH;
	    uint ygridlength = ydim/BLOCK_LENGTH;
	    uint zgridlength = zdim/BLOCK_LENGTH;

	    // compute block index
	    uint bx = block_idx%xgridlength;
	    uint tmpIdx = (block_idx - bx)/xgridlength;
	    uint by = tmpIdx%ygridlength;
	    uint bz = (tmpIdx-by)/ygridlength;

	    // copy global to shared memory
	    dim3 idx(tx+1,ty+1,tz+1);
			
            _sol[idx.x][idx.y][idx.z] = sol_in[base_addr + tIdx];
			
            float F = spd[base_addr + tIdx];
			
            bool isValid = mask[base_addr + tIdx];

	    uint new_base_addr, new_tIdx;

	    // 1-neighborhood values
	    if(tx == 0) 
	    {
		if(bx == 0) // end of the grid
		{	
	  	    new_tIdx = tIdx;
		    new_base_addr = base_addr;
		}
		else
		{
		    new_tIdx = tIdx + BLOCK_LENGTH-1;
		    new_base_addr = (block_idx - 1)*blocksize;	
		}
		_sol[tx][idx.y][idx.z] = sol_in[new_base_addr + new_tIdx];	
	    }

	    if(tx == BLOCK_LENGTH-1)
	    {
		if(bx == xgridlength-1) // end of the grid
		{
		    new_tIdx = tIdx;
		    new_base_addr = base_addr;
		}
		else
		{
		    new_tIdx = tIdx - (BLOCK_LENGTH-1);
		    new_base_addr = (block_idx + 1)*blocksize;	
		}
		_sol[tx+2][idx.y][idx.z] = sol_in[new_base_addr + new_tIdx];	
	    }

	    if(ty == 0)
	    {
		if(by == 0)
		{
		    new_tIdx = tIdx;
		    new_base_addr = base_addr;
		}
		else
		{
		    new_tIdx = tIdx + (BLOCK_LENGTH-1)*BLOCK_LENGTH;
		    new_base_addr = (block_idx - xgridlength)*blocksize;
		}
		_sol[idx.x][ty][idx.z] = sol_in[new_base_addr + new_tIdx];
	    }

	    if(ty == BLOCK_LENGTH-1) 
	    {
		if(by == ygridlength-1) 
		{
		    new_tIdx = tIdx;
		    new_base_addr = base_addr;
		}
		else
		{
		    new_tIdx = tIdx - (BLOCK_LENGTH-1)*BLOCK_LENGTH;
		    new_base_addr = (block_idx + xgridlength)*blocksize;
		}
		_sol[idx.x][ty+2][idx.z] = sol_in[new_base_addr + new_tIdx];
	    }

	    if(tz == 0)
	    {
		if(bz == 0)
		{
		    new_tIdx = tIdx;
		    new_base_addr = base_addr;
		}
		else
		{
		    new_tIdx = tIdx + (BLOCK_LENGTH-1)*BLOCK_LENGTH*BLOCK_LENGTH;
		    new_base_addr = (block_idx - xgridlength*ygridlength)*blocksize;
		}
		    _sol[idx.x][idx.y][tz] = sol_in[new_base_addr + new_tIdx];
	    }

	    if(tz == BLOCK_LENGTH-1)
	    {
	        if(bz == zgridlength-1) // end of the grid
		{
		    new_tIdx = tIdx;
		    new_base_addr = base_addr;
		}
		else
		{
		    new_tIdx = tIdx - (BLOCK_LENGTH-1)*BLOCK_LENGTH*BLOCK_LENGTH;
		    new_base_addr = (block_idx + xgridlength*ygridlength)*blocksize;
		}
		_sol[idx.x][idx.y][tz+2] = sol_in[new_base_addr + new_tIdx];
	    }

	    __syncthreads();


	    float a, b, c, oldT, newT;

	    // compute new value
	    oldT = newT = _sol[idx.x][idx.y][idx.z];

	    if(isValid)
	    {
		a = min(_sol[tx][idx.y][idx.z],_sol[tx+2][idx.y][idx.z]);
		b = min(_sol[idx.x][ty][idx.z],_sol[idx.x][ty+2][idx.z]);
		c = min(_sol[idx.x][idx.y][tz],_sol[idx.x][idx.y][tz+2]);

		float tmp = get_time_eikonal(a, b, c, dh, F);
		newT = min(tmp,oldT);

		sol_out[base_addr + tIdx] = newT;
	    }

	    // write back to global memory
	    float residue = oldT - newT;
	    con[base_addr + tIdx] = (residue < EPS) ? true : false;	
	}
    }
}

__device__ float get_time_eikonal(float a, float b, float c, float h, float s)
{
    float T_ijk, tmp;

    // a > b > c
    if(a < b) { tmp = a; a = b; b = tmp; }
    if(b < c) { tmp = b; b = c; c = tmp; }
    if(a < b) { tmp = a; a = b; b = tmp; }

    T_ijk = INF;

    if(c < INF)
    {
 	T_ijk = c + h*s;
		
	if(T_ijk > b) 
	{	
	    tmp = ((b + c) + sqrtf(2.0f*s*s*h*h - (b - c)*(b - c)))*0.5f;
	
	    if(tmp > b) T_ijk = tmp; 

	    if(T_ijk > a)	
	    {				
                tmp = (a + b + c)/3.0f + sqrtf(2.0f*(a*(b - a) + b*(c - b) + c*(a - c)) + 3.0f*s*s*h*h) / 3.0f;

		if(tmp > a) T_ijk = tmp;
	    }
	}
    }

    return T_ijk;
} 
