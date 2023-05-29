# include "FSM.cuh"

void Eikonal_fsm::parameters()
{
    padb = 1;

    nSweeps = 8;
    meshDim = 3;

    nxx = nx + 2*padb;
    nyy = ny + 2*padb;
    nzz = nz + 2*padb;

    threadsPerBlock = 256;

	totalLevels = (nxx - 1) + (nyy - 1) + (nzz - 1);

    title = "Eikonal solver for acoustic isotropic media\n\nSolving eikonal equation with the \033[32mNoble, Gesret and Belayouni (2014)\033[0;0m formulation\n";    
}

void Eikonal_fsm::components() 
{ 
	dh2i = 1.0f / (dh*dh);
	dh4i = dh2i * dh2i;
	dsum = 3 * dh2i;    

    int sgnv[nSweeps][meshDim] = {{1,1,1}, {0,1,1}, {1,1,0}, {0,1,0}, {1,0,1}, {0,0,1}, {1,0,0}, {0,0,0}};
    int sgnt[nSweeps][meshDim] = {{1,1,1}, {-1,1,1}, {1,1,-1}, {-1,1,-1}, {1,-1,1}, {-1,-1,1}, {1,-1,-1}, {-1,-1,-1}};

	int * h_sgnv = new int [nSweeps * meshDim]();
	int * h_sgnt = new int [nSweeps * meshDim](); 

	for (int index = 0; index < nSweeps * meshDim; index++)
	{
		int j = index / nSweeps;
		int i = index % nSweeps;				

		h_sgnv[i + j * nSweeps] = sgnv[i][j];
		h_sgnt[i + j * nSweeps] = sgnt[i][j];
	}

	cudaMalloc((void**)&(d_T), volsize*sizeof(float));
	cudaMalloc((void**)&(d_S), volsize*sizeof(float));

	cudaMalloc((void**)&(d_sgnv), nSweeps*meshDim*sizeof(int));
	cudaMalloc((void**)&(d_sgnt), nSweeps*meshDim*sizeof(int));

	cudaMemcpy(d_sgnv, h_sgnv, nSweeps*meshDim*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_sgnt, h_sgnt, nSweeps*meshDim*sizeof(int), cudaMemcpyHostToDevice);

	cudaMemcpy(d_S, S, volsize*sizeof(float), cudaMemcpyHostToDevice);

    delete[] h_sgnt;
    delete[] h_sgnv;
}

void Eikonal_fsm::initial_setup()
{
    int sidx = (int)(geometry->shots.x[shot_id] / dh) + padb;
    int sidy = (int)(geometry->shots.y[shot_id] / dh) + padb;
    int sidz = (int)(geometry->shots.z[shot_id] / dh) + padb;

    int sId = sidz + sidx*nzz + sidy*nxx*nzz;

    for (int index = 0; index < volsize; index++)
        T[index] = 1e6f;

    // Neighboring source points initialization with analitical traveltime

    T[sId] = S[sId] * sqrtf(powf((sidx-padb)*dh - geometry->shots.x[shot_id], 2.0f) + powf((sidy-padb)*dh - geometry->shots.y[shot_id], 2.0f) + powf((sidz-padb)*dh - geometry->shots.z[shot_id], 2.0f));
    
    T[sId + 1] = S[sId] * sqrtf(powf((sidx-padb)*dh - geometry->shots.x[shot_id], 2.0f) + powf((sidy-padb)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)+1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - 1] = S[sId] * sqrtf(powf((sidx-padb)*dh - geometry->shots.x[shot_id], 2.0f) + powf((sidy-padb)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)-1)*dh - geometry->shots.z[shot_id], 2.0f));

    T[sId + nzz] = S[sId] * sqrtf(powf(((sidx-padb)+1)*dh - geometry->shots.x[shot_id], 2.0f) + powf((sidy-padb)*dh - geometry->shots.y[shot_id], 2.0f) + powf((sidz-padb)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - nzz] = S[sId] * sqrtf(powf(((sidx-padb)-1)*dh - geometry->shots.x[shot_id], 2.0f) + powf((sidy-padb)*dh - geometry->shots.y[shot_id], 2.0f) + powf((sidz-padb)*dh - geometry->shots.z[shot_id], 2.0f));
    
    T[sId + nxx*nzz] = S[sId] * sqrtf(powf((sidx-padb)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)+1)*dh - geometry->shots.y[shot_id], 2.0f) + powf((sidz-padb)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - nxx*nzz] = S[sId] * sqrtf(powf((sidx-padb)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)-1)*dh - geometry->shots.y[shot_id], 2.0f) + powf((sidz-padb)*dh - geometry->shots.z[shot_id], 2.0f));
    
    T[sId + 1 + nzz] = S[sId] * sqrtf(powf(((sidx-padb)+1)*dh - geometry->shots.x[shot_id], 2.0f) + powf((sidy-padb)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)+1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId + 1 - nzz] = S[sId] * sqrtf(powf(((sidx-padb)+1)*dh - geometry->shots.x[shot_id], 2.0f) + powf((sidy-padb)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)-1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - 1 + nzz] = S[sId] * sqrtf(powf(((sidx-padb)-1)*dh - geometry->shots.x[shot_id], 2.0f) + powf((sidy-padb)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)+1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - 1 - nzz] = S[sId] * sqrtf(powf(((sidx-padb)-1)*dh - geometry->shots.x[shot_id], 2.0f) + powf((sidy-padb)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)-1)*dh - geometry->shots.z[shot_id], 2.0f));
    
    T[sId + 1 + nxx*nzz] = S[sId] * sqrtf(powf((sidx-padb)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)+1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)+1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId + 1 - nxx*nzz] = S[sId] * sqrtf(powf((sidx-padb)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)-1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)+1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - 1 + nxx*nzz] = S[sId] * sqrtf(powf((sidx-padb)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)+1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)-1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - 1 - nxx*nzz] = S[sId] * sqrtf(powf((sidx-padb)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)-1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)-1)*dh - geometry->shots.z[shot_id], 2.0f));
    
    T[sId + nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)+1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)+1)*dh - geometry->shots.y[shot_id], 2.0f) + powf((sidz-padb)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId + nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)+1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)-1)*dh - geometry->shots.y[shot_id], 2.0f) + powf((sidz-padb)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)-1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)+1)*dh - geometry->shots.y[shot_id], 2.0f) + powf((sidz-padb)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)-1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)-1)*dh - geometry->shots.y[shot_id], 2.0f) + powf((sidz-padb)*dh - geometry->shots.z[shot_id], 2.0f));
    
    T[sId + 1 + nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)+1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)+1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)+1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId + 1 - nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)-1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)+1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)+1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId + 1 + nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)+1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)-1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)+1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId + 1 - nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)-1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)-1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)+1)*dh - geometry->shots.z[shot_id], 2.0f));

    T[sId - 1 + nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)+1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)+1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)-1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - 1 - nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)-1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)+1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)-1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - 1 + nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)+1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)-1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)-1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - 1 - nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)-1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)-1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)-1)*dh - geometry->shots.z[shot_id], 2.0f));

	cudaMemcpy(d_T, T, nxx*nyy*nzz*sizeof(float), cudaMemcpyHostToDevice);
}

void Eikonal_fsm::expansion()
{
    for (int z = padb; z < nzz - padb; z++)
    {
        for (int y = padb; y < nyy - padb; y++)
        {
            for (int x = padb; x < nxx - padb; x++)
            {
                S[z + x*nzz + y*nxx*nzz] = 1.0f / V[(z - padb) + (x - padb)*nz + (y - padb)*nx*nz];
            }
        }
    }

    for (int z = 0; z < padb; z++)
    {
        for (int y = padb; y < nyy - padb; y++)
        {
            for (int x = padb; x < nxx - padb; x++)
            {
                S[z + x*nzz + y*nxx*nzz] = 1.0f / V[0 + (x - padb)*nz + (y - padb)*nx*nz];
                S[(nzz - z - 1) + x*nzz + y*nxx*nzz] = 1.0f / V[(nz - 1) + (x - padb)*nz + (y - padb)*nx*nz];
            }
        }
    }

    for (int x = 0; x < padb; x++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int y = padb; y < nyy - padb; y++)
            {
                S[z + x*nzz + y*nxx*nzz] = S[z + padb*nzz + y*nxx*nzz];
                S[z + (nxx - x - 1)*nzz + y*nxx*nzz] = S[z + (nxx - padb - 1)*nzz + y*nxx*nzz];
            }
        }
    }

    for (int y = 0; y < padb; y++)
    {
        for (int z = 0; z < nzz; z++)
        {
            for (int x = 0; x < nxx; x++)
            {
                S[z + x*nzz + y*nxx*nzz] = S[z + x*nzz + padb*nxx*nzz];
                S[z + x*nzz + (nyy - y - 1)*nxx*nzz] = S[z + x*nzz + (nyy - padb - 1)*nxx*nzz];
            }
        }
    }
}

void Eikonal_fsm::reduction()
{
    for (int index = 0; index < nPoints; index++)
    {
        int y = (int) (index / (nx*nz));         
        int x = (int) (index - y*nx*nz) / nz;    
        int z = (int) (index - x*nz - y*nx*nz);  

        wavefield_output[z + x*nz + y*nx*nz] = T[(z + padb) + (x + padb)*nzz + (y + padb)*nxx*nzz];
    }
}

void Eikonal_fsm::forward_solver()
{    
    for (int sweep = 0; sweep < 8; sweep++)
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
			
			fast_sweeping_kernel<<<gs,bs>>>(d_S,d_T,d_sgnt,d_sgnv,sweep,nSweeps,level,xs,ys,xSweepOff,ySweepOff,zSweepOff,nxx,nyy,nzz,dh,dh2i,dh4i,dsum);
			
            cudaDeviceSynchronize();
		}
	}

    cudaMemcpy(T, d_T, volsize*sizeof(float), cudaMemcpyDeviceToHost);
}

int Eikonal_fsm::iDivUp(int a, int b) { return ( (a % b) != 0 ) ? (a / b + 1) : (a / b); }

void Eikonal_fsm::free_space()
{
    cudaFree(d_T);
    cudaFree(d_S);

    cudaFree(d_sgnt);
    cudaFree(d_sgnv);

    delete[] T;
    delete[] S;
}

__global__ void fast_sweeping_kernel(float * S, float * T, int * sgnt, int * sgnv, int sc, int nSweeps, 
                                     int level, int xOffset, int yOffset, int xSweepOffset, int ySweepOffset, 
                                     int zSweepOffset, int nxx, int nyy, int nzz, float dh, float dh2i, 
                                     float dh4i, float dsum)
{
	int x = (blockIdx.x * blockDim.x + threadIdx.x) + xOffset;
	int y = (blockIdx.y * blockDim.y + threadIdx.y) + yOffset;

    float ta, tb, tc, t1, t2, t3, Sref;
    float t1D1, t1D2, t1D3, t1D, t2D1, t2D2, t2D3, t2D, t3D;

	if ((x > 0) && (x <= nxx) && (y > 0) && (y <= nyy)) 
	{
		int z = level - (x + y);
		
		if ((z > 0) && (z <= nzz))	
		{
			int i = abs(z - zSweepOffset);
			int j = abs(x - xSweepOffset);
			int k = abs(y - ySweepOffset);

			if ((i > 0) && (i < nzz-1) && (j > 0) && (j < nxx-1) && (k > 0) && (k < nyy-1))
			{		
				int i1 = i - sgnv[sc + 0*nSweeps];
				int j1 = j - sgnv[sc + 1*nSweeps];
				int k1 = k - sgnv[sc + 2*nSweeps];

				int ijk = i + j*nzz + k*nxx*nzz;
				
				float tv = T[(i - sgnt[sc + 0*nSweeps]) + j*nzz + k*nxx*nzz];
				float te = T[i + (j - sgnt[sc + 1*nSweeps])*nzz + k*nxx*nzz];
				float tn = T[i + j*nzz + (k - sgnt[sc + 2*nSweeps])*nxx*nzz];

				float tev = T[(i - sgnt[sc + 0*nSweeps]) + (j - sgnt[sc + 1*nSweeps])*nzz + k*nxx*nzz];
				float ten = T[i + (j - sgnt[sc + 1*nSweeps])*nzz + (k - sgnt[sc + 2*nSweeps])*nxx*nzz];
				float tnv = T[(i - sgnt[sc + 0*nSweeps]) + j*nzz + (k - sgnt[sc + 2*nSweeps])*nxx*nzz];
				
				float tnve = T[(i - sgnt[sc + 0*nSweeps]) + (j - sgnt[sc + 1*nSweeps])*nzz + (k - sgnt[sc + 2*nSweeps])*nxx*nzz];

				t1D1 = tv + dh * min(S[i1 + max(j-1,1)*nzz   + max(k-1,1)*nxx*nzz], 
								 min(S[i1 + max(j-1,1)*nzz   + min(k,nyy-1)*nxx*nzz], 
								 min(S[i1 + min(j,nxx-1)*nzz + max(k-1,1)*nxx*nzz],
									 S[i1 + min(j,nxx-1)*nzz + min(k,nyy-1)*nxx*nzz])));                                     

				t1D2 = te + dh * min(S[max(i-1,1)   + j1*nzz + max(k-1,1)*nxx*nzz], 
								 min(S[min(i,nzz-1) + j1*nzz + max(k-1,1)*nxx*nzz],
								 min(S[max(i-1,1)   + j1*nzz + min(k,nyy-1)*nxx*nzz], 
									 S[min(i,nzz-1) + j1*nzz + min(k,nyy-1)*nxx*nzz])));                    

				t1D3 = tn + dh * min(S[max(i-1,1)   + max(j-1,1)*nzz   + k1*nxx*nzz], 
								 min(S[max(i-1,1)   + min(j,nxx-1)*nzz + k1*nxx*nzz],
								 min(S[min(i,nzz-1) + max(j-1,1)*nzz   + k1*nxx*nzz], 
									 S[min(i,nzz-1) + min(j,nxx-1)*nzz + k1*nxx*nzz])));

				t1D = min(t1D1, min(t1D2, t1D3));

				//------------------- 2D operators - 4 points operator ---------------------------------------------------------------------------------------------------
				t2D1 = 1e6f; t2D2 = 1e6f; t2D3 = 1e6f;

				// XZ plane ----------------------------------------------------------------------------------------------------------------------------------------------
				Sref = min(S[i1 + j1*nzz + max(k-1,1)*nxx*nzz], S[i1 + j1*nzz + min(k, nyy-1)*nxx*nzz]);
				
				if ((tv < te + dh*Sref) && (te < tv + dh*Sref))
				{
					ta = tev + te - tv;
					tb = tev - te + tv;

					t2D1 = ((tb*dh2i + ta*dh2i) + sqrtf(4.0f*Sref*Sref*(dh2i + dh2i) - dh2i*dh2i*(ta - tb)*(ta - tb))) / (dh2i + dh2i);
				}

				// YZ plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
				Sref = min(S[i1 + max(j-1,1)*nzz + k1*nxx*nzz], S[i1 + min(j,nxx-1)*nzz + k1*nxx*nzz]);

				if((tv < tn + dh*Sref) && (tn < tv + dh*Sref))
				{
					ta = tv - tn + tnv;
					tb = tn - tv + tnv;
					
					t2D2 = ((ta*dh2i + tb*dh2i) + sqrtf(4.0f*Sref*Sref*(dh2i + dh2i) - dh2i*dh2i*(ta - tb)*(ta - tb))) / (dh2i + dh2i); 
				}

				// XY plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
				Sref = min(S[max(i-1,1) + j1*nzz + k1*nxx*nzz], S[min(i,nzz-1) + j1*nzz + k1*nxx*nzz]);

				if((te < tn + dh*Sref) && (tn < te + dh*Sref))
				{
					ta = te - tn + ten;
					tb = tn - te + ten;

					t2D3 = ((ta*dh2i + tb*dh2i) + sqrtf(4.0f*Sref*Sref*(dh2i + dh2i) - dh2i*dh2i*(ta - tb)*(ta - tb))) / (dh2i + dh2i);
				}

				t2D = min(t2D1, min(t2D2, t2D3));

				//------------------- 3D operators ---------------------------------------------------------------------------------------------------
				t3D = 1e6f;

				Sref = S[i1 + j1*nzz + k1*nxx*nzz];

				ta = te - 0.5f*tn + 0.5f*ten - 0.5f*tv + 0.5f*tev - tnv + tnve;
				tb = tv - 0.5f*tn + 0.5f*tnv - 0.5f*te + 0.5f*tev - ten + tnve;
				tc = tn - 0.5f*te + 0.5f*ten - 0.5f*tv + 0.5f*tnv - tev + tnve;

				if (min(t1D,t2D) > max(tv, max(te, tn)))
				{
					t2 = 9.0f*Sref*Sref*dsum;
					
					t3 = dh4i*(ta - tb)*(ta - tb) + dh4i*(tb - tc)*(tb - tc) + dh4i*(ta - tc)*(ta - tc);
					
					if (t2 >= t3)
					{
						t1 = tb*dh2i + ta*dh2i + tc*dh2i;        
						
						t3D = (t1 + sqrtf(t2 - t3)) / dsum;
					}
				}

				T[ijk] = min(T[ijk], min(t1D, min(t2D, t3D)));
            }
        }
    }
}
