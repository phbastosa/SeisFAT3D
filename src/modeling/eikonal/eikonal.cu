# include "eikonal.cuh"

void Eikonal::set_parameters(std::string file)
{
    Modeling::set_parameters(file);

    std::string vp_model_file = catch_parameter("vp_model_file", file);

    // export_receiver_output = 
    // export_wavefield_output = 

    // receiver_output_folder = 
    // wavefield_output_folder = 

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

    wavefield_output_samples = nPoints;
    receiver_output_samples = geometry->nodes.total;

    receiver_output = new float[receiver_output_samples]();
    wavefield_output = new float[wavefield_output_samples]();
}

void Eikonal::initial_setup()
{
	uint idx = 0;
	uint blk_idx = 0;
	uint list_idx = 0;
    uint nActiveBlock = 0;

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
}

void Eikonal::forward_solver()
{







}

void Eikonal::build_outputs()
{
    // receiver_output_file = receiver_output_folder + " ";
    // wavefield_output_file = wavefield_output_file + " ";

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
}