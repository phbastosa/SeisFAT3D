# include "eikonal.cuh"

void Eikonal::set_parameters(std::string file)
{
    Modeling::set_parameters(file);

    title = "Eikonal solver for acoustic isotropic media\n";

    std::string vp_model_file = catch_parameter("vp_model_file", file);

    pdx = (BLOCK_LENGTH - nx % BLOCK_LENGTH) % BLOCK_LENGTH;
    pdy = (BLOCK_LENGTH - ny % BLOCK_LENGTH) % BLOCK_LENGTH;
    pdz = (BLOCK_LENGTH - nz % BLOCK_LENGTH) % BLOCK_LENGTH;

    nxx = nx + pdx;    
    nyy = ny + pdy;    
    nzz = nz + pdz;    
    
	volsize = nxx * nyy * nzz;

	blksize = BLOCK_LENGTH * BLOCK_LENGTH * BLOCK_LENGTH;

	nbx = nxx / BLOCK_LENGTH;
	nby = nyy / BLOCK_LENGTH;
	nbz = nzz / BLOCK_LENGTH;
	
    nblk = nbx * nby * nbz;

    S = new float[volsize]();
    vp = new float[nPoints]();

    import_binary_float(vp_model_file, vp, nPoints);

    expand_model();

    delete[] vp;
}

void Eikonal::expand_model()
{
    for (int z = 0; z < nz; z++)
    {
        for (int y = 0; y < ny; y++)
        {
            for (int x = 0; x < nx; x++)
            {
                S[z + x*nzz + y*nxx*nzz] = 1.0f / vp[z + x*nz + y*nx*nz];
            }
        }
    }

    for (int z = 0; z < pdz; z++)
    {
        for (int y = 0; y < nyy - pdy; y++)
        {
            for (int x = 0; x < nxx - pdx; x++)
            {
                S[(nzz - z - 1) + x*nzz + y*nxx*nzz] = 1.0f / vp[(nz - 1) + x*nz + y*nx*nz];
            }
        }
    }

    for (int x = 0; x < pdx; x++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int y = 0; y < nyy - pdy; y++)
            {
                S[z + (nxx - x - 1)*nzz + y*nxx*nzz] = S[z + (nxx - pdx - 1)*nzz + y*nxx*nzz];
            }
        }
    }

    for (int y = 0; y < pdy; y++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int x = 0; x < nxx; x++)
            {
                S[z + x*nzz + (nyy - y - 1)*nxx*nzz] = S[z + x*nzz + (nyy - pdy - 1)*nxx*nzz];
            }
        }
    }
}

void Eikonal::set_components()
{
    T = new float[nPoints]();

    h_mask = new bool[volsize]();
    h_slow = new float[volsize]();
	h_time = new float[volsize]();
	
    h_list = new uint[nblk]();      
	h_listed = new bool[nblk]();    
	h_listVol = new bool[nblk]();   

	cudaMalloc((void**)&(d_con), volsize*sizeof(bool));  
	cudaMalloc((void**)&(d_mask), volsize*sizeof(bool));
	cudaMalloc((void**)&(d_slow), volsize*sizeof(float));
	cudaMalloc((void**)&(d_time), volsize*sizeof(float));
	cudaMalloc((void**)&(t_time), volsize*sizeof(float)); 
	cudaMalloc((void**)&(d_list), nblk*sizeof(uint));
	cudaMalloc((void**)&(d_listVol), nblk*sizeof(bool));

    wavefield_output_samples = nPoints;
    receiver_output_samples = geometry->nodes.total;

    receiver_output = new float[receiver_output_samples]();
    wavefield_output = new float[wavefield_output_samples]();
    
    get_RAM_usage();
    get_GPU_usage();
}

void Eikonal::initial_setup()
{
	uint idx = 0;
	uint blk_idx = 0;
	uint list_idx = 0;
    nActiveBlock = 0;

    float sx = geometry->shots.x[shot_id];
    float sy = geometry->shots.y[shot_id];
    float sz = geometry->shots.z[shot_id];

    int sidx = (int)(sx / dh);
    int sidy = (int)(sy / dh);
    int sidz = (int)(sz / dh);

    int sid = sidz + sidx*nzz + sidy*nxx*nzz;

    float t0 = S[sid] * sqrtf(powf((float)(sidx*dh) - sx, 2.0f) +
                              powf((float)(sidy*dh) - sy, 2.0f) +
                              powf((float)(sidz*dh) - sz, 2.0f));

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

void Eikonal::forward_solver()
{
    FIM_solver();
}

void Eikonal::build_outputs()
{
    get_travelTimes();
    get_firstArrivals();
}

void Eikonal::get_travelTimes()
{
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
							if ((z < nz) && (x < nx) && (y < ny))
                            {
                                T[z + x*nz + y*nx*nz] = h_time[idx];
                                wavefield_output[z + x*nz + y*nx*nz] = h_time[idx];
                            } 
			
							++idx;
						}
					}
				}
			}
		}
	}

    wavefield_output_file = wavefield_output_folder + "travel_times_" + std::to_string(nz) + "x" + std::to_string(nx) + "x" + std::to_string(ny) + "_shot_" + std::to_string(shot_id+1) + ".bin";
}

void Eikonal::get_firstArrivals()
{
    for (int r = 0; r < total_nodes; r++)
    {
        float x = geometry->nodes.x[r];
        float y = geometry->nodes.y[r];
        float z = geometry->nodes.z[r];

        float x0 = floorf(x / dh) * dh;
        float y0 = floorf(y / dh) * dh;
        float z0 = floorf(z / dh) * dh;

        float x1 = floorf(x / dh) * dh + dh;
        float y1 = floorf(y / dh) * dh + dh;
        float z1 = floorf(z / dh) * dh + dh;

        int id = ((int)(z / dh)) + ((int)(x / dh))*nz + ((int)(y / dh))*nx*nz;

        float c000 = T[id];
        float c001 = T[id + 1];
        float c100 = T[id + nz]; 
        float c101 = T[id + 1 + nz]; 
        float c010 = T[id + nx*nz]; 
        float c011 = T[id + 1 + nx*nz]; 
        float c110 = T[id + nz + nx*nz]; 
        float c111 = T[id + 1 + nz + nx*nz];

        float xd = (x - x0) / (x1 - x0);
        float yd = (y - y0) / (y1 - y0);
        float zd = (z - z0) / (z1 - z0);

        float c00 = c000*(1 - xd) + c100*xd;    
        float c01 = c001*(1 - xd) + c101*xd;    
        float c10 = c010*(1 - xd) + c110*xd;    
        float c11 = c011*(1 - xd) + c111*xd;    

        float c0 = c00*(1 - yd) + c10*yd;
        float c1 = c01*(1 - yd) + c11*yd;

        receiver_output[r] = c0*(1 - zd) + c1*zd;
    }

    receiver_output_file = receiver_output_folder + "first_arrivals_" + std::to_string(geometry->nodes.total) + "_shot_" + std::to_string(shot_id+1) + ".bin";
}

void Eikonal::FIM_solver()
{
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
