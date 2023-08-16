# include "fast_sweeping_method.cuh"

void Fast_Sweeping_Method::set_parameters()
{
    general_modeling_parameters();

    set_acquisition_geometry();

    set_velocity_model();
    
    set_boundaries();
    set_model_boundaries();
    
    set_slowness_model();
    set_outputs();

    set_modeling_volumes();
}

void Fast_Sweeping_Method::set_boundaries()
{
    nbxl = 1; nbxr = 1;
    nbyl = 1; nbyr = 1;
    nbzu = 1; nbzd = 1;
}

void Fast_Sweeping_Method::set_modeling_volumes()
{
    modeling_method = std::string("fsm");

    nSweeps = 8;
    meshDim = 3;

    dz2i = 1.0f / (dz*dz);
    dx2i = 1.0f / (dx*dx);
    dy2i = 1.0f / (dy*dy);

    dz2dx2 = dz2i * dx2i;
    dz2dy2 = dz2i * dy2i;
    dx2dy2 = dx2i * dy2i;

    dsum = dz2i + dx2i + dy2i;

	totalLevels = (nxx - 1) + (nyy - 1) + (nzz - 1);

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

    T = new float[volsize](); 

	cudaMalloc((void**)&(d_T), volsize*sizeof(float));
	cudaMalloc((void**)&(d_S), volsize*sizeof(float));

	cudaMalloc((void**)&(d_sgnv), nSweeps*meshDim*sizeof(int));
	cudaMalloc((void**)&(d_sgnt), nSweeps*meshDim*sizeof(int));

	cudaMemcpy(d_sgnv, h_sgnv, nSweeps*meshDim*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sgnt, h_sgnt, nSweeps*meshDim*sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy(d_S, S, volsize*sizeof(float), cudaMemcpyHostToDevice);

    delete[] h_sgnt;
    delete[] h_sgnv;

    std::vector<std::vector<int>>().swap(sgnv);
    std::vector<std::vector<int>>().swap(sgnt);
}

void Fast_Sweeping_Method::info_message()
{
    general_modeling_message();

    std::cout<<"[2] - Fast Sweeping Method (Detrixhe et al., 2013; Noble et al., 2014)\n\n"; 
}

void Fast_Sweeping_Method::initial_setup()
{
    int sidx = (int)(geometry->shots.x[shot_id] / dx) + nbxl;
    int sidy = (int)(geometry->shots.y[shot_id] / dy) + nbyl;
    int sidz = (int)(geometry->shots.z[shot_id] / dz) + nbzu;

    source_id = sidz + sidx*nzz + sidy*nxx*nzz;

    for (int index = 0; index < volsize; index++) T[index] = 1e6f;

    T[source_id] = S[source_id] * sqrtf(powf((sidx-nbxl)*dx - geometry->shots.x[shot_id], 2.0f) + powf((sidy-nbyl)*dy - geometry->shots.y[shot_id], 2.0f) + powf((sidz-nbzu)*dz - geometry->shots.z[shot_id], 2.0f));

    T[source_id + 1] = S[source_id] * sqrtf(powf((sidx-nbxl)*dx - geometry->shots.x[shot_id], 2.0f) + powf((sidy-nbyl)*dy - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-nbzu)+1)*dz - geometry->shots.z[shot_id], 2.0f));
    T[source_id - 1] = S[source_id] * sqrtf(powf((sidx-nbxl)*dx - geometry->shots.x[shot_id], 2.0f) + powf((sidy-nbyl)*dy - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-nbzu)-1)*dz - geometry->shots.z[shot_id], 2.0f));

    T[source_id + nzz] = S[source_id] * sqrtf(powf(((sidx-nbxl)+1)*dx - geometry->shots.x[shot_id], 2.0f) + powf((sidy-nbyl)*dy - geometry->shots.y[shot_id], 2.0f) + powf((sidz-nbzu)*dz - geometry->shots.z[shot_id], 2.0f));
    T[source_id - nzz] = S[source_id] * sqrtf(powf(((sidx-nbxl)-1)*dx - geometry->shots.x[shot_id], 2.0f) + powf((sidy-nbyl)*dy - geometry->shots.y[shot_id], 2.0f) + powf((sidz-nbzu)*dz - geometry->shots.z[shot_id], 2.0f));
    
    T[source_id + nxx*nzz] = S[source_id] * sqrtf(powf((sidx-nbxl)*dx - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-nbyl)+1)*dy - geometry->shots.y[shot_id], 2.0f) + powf((sidz-nbzu)*dz - geometry->shots.z[shot_id], 2.0f));
    T[source_id - nxx*nzz] = S[source_id] * sqrtf(powf((sidx-nbxl)*dx - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-nbyl)-1)*dy - geometry->shots.y[shot_id], 2.0f) + powf((sidz-nbzu)*dz - geometry->shots.z[shot_id], 2.0f));
    
    T[source_id + 1 + nzz] = S[source_id] * sqrtf(powf(((sidx-nbxl)+1)*dx - geometry->shots.x[shot_id], 2.0f) + powf((sidy-nbyl)*dy - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-nbzu)+1)*dz - geometry->shots.z[shot_id], 2.0f));
    T[source_id + 1 - nzz] = S[source_id] * sqrtf(powf(((sidx-nbxl)+1)*dx - geometry->shots.x[shot_id], 2.0f) + powf((sidy-nbyl)*dy - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-nbzu)-1)*dz - geometry->shots.z[shot_id], 2.0f));
    T[source_id - 1 + nzz] = S[source_id] * sqrtf(powf(((sidx-nbxl)-1)*dx - geometry->shots.x[shot_id], 2.0f) + powf((sidy-nbyl)*dy - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-nbzu)+1)*dz - geometry->shots.z[shot_id], 2.0f));
    T[source_id - 1 - nzz] = S[source_id] * sqrtf(powf(((sidx-nbxl)-1)*dx - geometry->shots.x[shot_id], 2.0f) + powf((sidy-nbyl)*dy - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-nbzu)-1)*dz - geometry->shots.z[shot_id], 2.0f));
    
    T[source_id + 1 + nxx*nzz] = S[source_id] * sqrtf(powf((sidx-nbxl)*dx - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-nbyl)+1)*dy - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-nbzu)+1)*dz - geometry->shots.z[shot_id], 2.0f));
    T[source_id + 1 - nxx*nzz] = S[source_id] * sqrtf(powf((sidx-nbxl)*dx - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-nbyl)-1)*dy - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-nbzu)+1)*dz - geometry->shots.z[shot_id], 2.0f));
    T[source_id - 1 + nxx*nzz] = S[source_id] * sqrtf(powf((sidx-nbxl)*dx - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-nbyl)+1)*dy - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-nbzu)-1)*dz - geometry->shots.z[shot_id], 2.0f));
    T[source_id - 1 - nxx*nzz] = S[source_id] * sqrtf(powf((sidx-nbxl)*dx - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-nbyl)-1)*dy - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-nbzu)-1)*dz - geometry->shots.z[shot_id], 2.0f));
    
    T[source_id + nzz + nxx*nzz] = S[source_id] * sqrtf(powf(((sidx-nbxl)+1)*dx - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-nbyl)+1)*dy - geometry->shots.y[shot_id], 2.0f) + powf((sidz-nbzu)*dz - geometry->shots.z[shot_id], 2.0f));
    T[source_id + nzz - nxx*nzz] = S[source_id] * sqrtf(powf(((sidx-nbxl)+1)*dx - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-nbyl)-1)*dy - geometry->shots.y[shot_id], 2.0f) + powf((sidz-nbzu)*dz - geometry->shots.z[shot_id], 2.0f));
    T[source_id - nzz + nxx*nzz] = S[source_id] * sqrtf(powf(((sidx-nbxl)-1)*dx - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-nbyl)+1)*dy - geometry->shots.y[shot_id], 2.0f) + powf((sidz-nbzu)*dz - geometry->shots.z[shot_id], 2.0f));
    T[source_id - nzz - nxx*nzz] = S[source_id] * sqrtf(powf(((sidx-nbxl)-1)*dx - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-nbyl)-1)*dy - geometry->shots.y[shot_id], 2.0f) + powf((sidz-nbzu)*dz - geometry->shots.z[shot_id], 2.0f));
    
    T[source_id + 1 + nzz + nxx*nzz] = S[source_id] * sqrtf(powf(((sidx-nbxl)+1)*dx - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-nbyl)+1)*dy - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-nbzu)+1)*dz - geometry->shots.z[shot_id], 2.0f));
    T[source_id + 1 - nzz + nxx*nzz] = S[source_id] * sqrtf(powf(((sidx-nbxl)-1)*dx - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-nbyl)+1)*dy - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-nbzu)+1)*dz - geometry->shots.z[shot_id], 2.0f));
    T[source_id + 1 + nzz - nxx*nzz] = S[source_id] * sqrtf(powf(((sidx-nbxl)+1)*dx - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-nbyl)-1)*dy - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-nbzu)+1)*dz - geometry->shots.z[shot_id], 2.0f));
    T[source_id + 1 - nzz - nxx*nzz] = S[source_id] * sqrtf(powf(((sidx-nbxl)-1)*dx - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-nbyl)-1)*dy - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-nbzu)+1)*dz - geometry->shots.z[shot_id], 2.0f));

    T[source_id - 1 + nzz + nxx*nzz] = S[source_id] * sqrtf(powf(((sidx-nbxl)+1)*dx - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-nbyl)+1)*dy - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-nbzu)-1)*dz - geometry->shots.z[shot_id], 2.0f));
    T[source_id - 1 - nzz + nxx*nzz] = S[source_id] * sqrtf(powf(((sidx-nbxl)-1)*dx - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-nbyl)+1)*dy - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-nbzu)-1)*dz - geometry->shots.z[shot_id], 2.0f));
    T[source_id - 1 + nzz - nxx*nzz] = S[source_id] * sqrtf(powf(((sidx-nbxl)+1)*dx - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-nbyl)-1)*dy - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-nbzu)-1)*dz - geometry->shots.z[shot_id], 2.0f));
    T[source_id - 1 - nzz - nxx*nzz] = S[source_id] * sqrtf(powf(((sidx-nbxl)-1)*dx - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-nbyl)-1)*dy - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-nbzu)-1)*dz - geometry->shots.z[shot_id], 2.0f));

    t0 = T[source_id];
}

void Fast_Sweeping_Method::forward_solver()
{
	cudaMemcpy(d_T, T, volsize*sizeof(float), cudaMemcpyHostToDevice);

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

			fast_sweeping_kernel<<<gs,bs>>>(d_S, d_T, d_sgnt, d_sgnv, sgni, sgnj, sgnk, level, xs, ys, 
                                            xSweepOff, ySweepOff, zSweepOff, nxx, nyy, nzz, dx, dy, dz, 
                                            dx2i, dy2i, dz2i, dz2dx2, dz2dy2, dx2dy2, dsum);
			
            cudaDeviceSynchronize();
		}
	}

    cudaMemcpy(T, d_T, volsize*sizeof(float), cudaMemcpyDeviceToHost);
}

int Fast_Sweeping_Method::iDivUp(int a, int b) 
{ 
    return ( (a % b) != 0 ) ? (a / b + 1) : (a / b); 
}

void Fast_Sweeping_Method::free_space()
{
    cudaFree(d_T);
    cudaFree(d_S);

    cudaFree(d_sgnt);
    cudaFree(d_sgnv);
}

__global__ void fast_sweeping_kernel(float * S, float * T, int * sgnt, int * sgnv, int sgni, int sgnj, int sgnk, 
                                     int level, int xOffset, int yOffset, int xSweepOffset, int ySweepOffset, int zSweepOffset, 
                                     int nxx, int nyy, int nzz, float dx, float dy, float dz, float dx2i, float dy2i, float dz2i, 
                                     float dz2dx2, float dz2dy2, float dx2dy2, float dsum)
{
	int x = (blockIdx.x * blockDim.x + threadIdx.x) + xOffset;
	int y = (blockIdx.y * blockDim.y + threadIdx.y) + yOffset;

    float ta, tb, tc, t1, t2, t3, Sref;
    float t1D1, t1D2, t1D3, t1D, t2D1, t2D2, t2D3, t2D, t3D;

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
				int i1 = i - sgnv[sgni];
				int j1 = j - sgnv[sgnj];
				int k1 = k - sgnv[sgnk];

				int ijk = i + j*nzz + k*nxx*nzz;
				
				float tv = T[(i - sgnt[sgni]) + j*nzz + k*nxx*nzz];
				float te = T[i + (j - sgnt[sgnj])*nzz + k*nxx*nzz];
				float tn = T[i + j*nzz + (k - sgnt[sgnk])*nxx*nzz];

				float tev = T[(i - sgnt[sgni]) + (j - sgnt[sgnj])*nzz + k*nxx*nzz];
				float ten = T[i + (j - sgnt[sgnj])*nzz + (k - sgnt[sgnk])*nxx*nzz];
				float tnv = T[(i - sgnt[sgni]) + j*nzz + (k - sgnt[sgnk])*nxx*nzz];
				
				float tnve = T[(i - sgnt[sgni]) + (j - sgnt[sgnj])*nzz + (k - sgnt[sgnk])*nxx*nzz];

				t1D1 = tv + dz * min(S[i1 + max(j-1,1)*nzz   + max(k-1,1)*nxx*nzz], 
								 min(S[i1 + max(j-1,1)*nzz   + min(k,nyy-1)*nxx*nzz], 
								 min(S[i1 + min(j,nxx-1)*nzz + max(k-1,1)*nxx*nzz],
									 S[i1 + min(j,nxx-1)*nzz + min(k,nyy-1)*nxx*nzz])));                                     

				t1D2 = te + dx * min(S[max(i-1,1)   + j1*nzz + max(k-1,1)*nxx*nzz], 
								 min(S[min(i,nzz-1) + j1*nzz + max(k-1,1)*nxx*nzz],
								 min(S[max(i-1,1)   + j1*nzz + min(k,nyy-1)*nxx*nzz], 
									 S[min(i,nzz-1) + j1*nzz + min(k,nyy-1)*nxx*nzz])));                    

				t1D3 = tn + dy * min(S[max(i-1,1)   + max(j-1,1)*nzz   + k1*nxx*nzz], 
								 min(S[max(i-1,1)   + min(j,nxx-1)*nzz + k1*nxx*nzz],
								 min(S[min(i,nzz-1) + max(j-1,1)*nzz   + k1*nxx*nzz], 
									 S[min(i,nzz-1) + min(j,nxx-1)*nzz + k1*nxx*nzz])));

				t1D = min(t1D1, min(t1D2, t1D3));

                //------------------- 2D operators - 4 points operator ---------------------------------------------------------------------------------------------------
                t2D1 = 1e6; t2D2 = 1e6; t2D3 = 1e6;

                // XZ plane ----------------------------------------------------------------------------------------------------------------------------------------------
                Sref = min(S[i1 + j1*nzz + max(k-1,1)*nxx*nzz], S[i1 + j1*nzz + min(k, nyy-1)*nxx*nzz]);
                
                if ((tv < te + dx*Sref) && (te < tv + dz*Sref))
                {
                    ta = tev + te - tv;
                    tb = tev - te + tv;

                    t2D1 = ((tb*dz2i + ta*dx2i) + sqrtf(4.0f*Sref*Sref*(dz2i + dx2i) - dz2i*dx2i*(ta - tb)*(ta - tb))) / (dz2i + dx2i);
                }

                // YZ plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
                Sref = min(S[i1 + max(j-1,1)*nzz + k1*nxx*nzz], S[i1 + min(j,nxx-1)*nzz + k1*nxx*nzz]);

                if((tv < tn + dy*Sref) && (tn < tv + dz*Sref))
                {
                    ta = tv - tn + tnv;
                    tb = tn - tv + tnv;
                    
                    t2D2 = ((ta*dz2i + tb*dy2i) + sqrtf(4.0f*Sref*Sref*(dz2i + dy2i) - dz2i*dy2i*(ta - tb)*(ta - tb))) / (dz2i + dy2i); 
                }

                // XY plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
                Sref = min(S[max(i-1,1) + j1*nzz + k1*nxx*nzz],S[min(i,nzz-1) + j1*nzz + k1*nxx*nzz]);

                if((te < tn + dy*Sref) && (tn < te + dx*Sref))
                {
                    ta = te - tn + ten;
                    tb = tn - te + ten;

                    t2D3 = ((ta*dx2i + tb*dy2i) + sqrtf(4.0f*Sref*Sref*(dx2i + dy2i) - dx2i*dy2i*(ta - tb)*(ta - tb))) / (dx2i + dy2i);
                }

                t2D = min(t2D1, min(t2D2, t2D3));

                //------------------- 3D operators - 8 point operator ---------------------------------------------------------------------------------------------------
                t3D = 1e6;

                Sref = S[i1 + j1*nzz + k1*nxx*nzz];

                ta = te - 0.5f*tn + 0.5f*ten - 0.5f*tv + 0.5f*tev - tnv + tnve;
                tb = tv - 0.5f*tn + 0.5f*tnv - 0.5f*te + 0.5f*tev - ten + tnve;
                tc = tn - 0.5f*te + 0.5f*ten - 0.5f*tv + 0.5f*tnv - tev + tnve;

                if (min(t1D, t2D) > max(tv, max(te, tn)))
                {
                    t2 = 9.0f*Sref*Sref*dsum; 
                    
                    t3 = dz2dx2*(ta - tb)*(ta - tb) + dz2dy2*(tb - tc)*(tb - tc) + dx2dy2*(ta - tc)*(ta - tc);
                    
                    if (t2 >= t3)
                    {
                        t1 = tb*dz2i + ta*dx2i + tc*dy2i;        
                        
                        t3D = (t1 + sqrtf(t2 - t3)) / dsum;
                    }
                }

				T[ijk] = min(T[ijk], min(t1D, min(t2D, t3D)));
            }
        }
    }
}
