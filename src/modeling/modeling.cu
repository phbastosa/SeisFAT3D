# include "modeling.cuh"

void Modeling::set_parameters()
{
    nx = std::stoi(catch_parameter("x_samples", parameters));
    ny = std::stoi(catch_parameter("y_samples", parameters));
    nz = std::stoi(catch_parameter("z_samples", parameters));

    dx = std::stof(catch_parameter("x_spacing", parameters));
    dy = std::stof(catch_parameter("y_spacing", parameters));
    dz = std::stof(catch_parameter("z_spacing", parameters));

    data_folder = catch_parameter("modeling_output_folder", parameters);

    nPoints = nx*ny*nz;

    geometry = new Geometry();
    geometry->parameters = parameters;
    geometry->set_parameters();

    max_spread = 0;
    for (int index = 0; index < geometry->nrel; index++)
    {   
        if (max_spread < geometry->spread[index])
            max_spread = geometry->spread[index]; 
    }

    seismogram = new float[max_spread]();

    nb = 3;
    
    nxx = nx + 2*nb;
    nyy = ny + 2*nb;
    nzz = nz + 2*nb;

    volsize = nxx*nyy*nzz;

    nThreads = 256;
    nBlocks = (int)((volsize + nThreads - 1) / nThreads);

    set_properties();    
    set_conditions();    
    set_eikonal();
}

void Modeling::set_properties()
{
    float * vp = new float[nPoints]();

    std::string vp_file = catch_parameter("vp_model_file", parameters);

    import_binary_float(vp_file, vp, nPoints);

    # pragma omp parallel for
    for (int index = 0; index < nPoints; index++)
        vp[index] = 1.0f / vp[index];

    S = new float[volsize]();

    expand_boundary(vp, S);

    delete[] vp;
}

void Modeling::set_eikonal()
{
    dz2i = 1.0f / (dz*dz);
    dx2i = 1.0f / (dx*dx);
    dy2i = 1.0f / (dy*dy);

    dz2dx2 = dz2i * dx2i;
    dz2dy2 = dz2i * dy2i;
    dx2dy2 = dx2i * dy2i;

    dsum = dz2i + dx2i + dy2i;

    total_levels = (nxx - 1) + (nyy - 1) + (nzz - 1);

    std::vector<std::vector<int>> sgnv = {{1,1,1}, {0,1,1}, {1,1,0}, {0,1,0}, {1,0,1}, {0,0,1}, {1,0,0}, {0,0,0}};
    std::vector<std::vector<int>> sgnt = {{1,1,1}, {-1,1,1}, {1,1,-1}, {-1,1,-1}, {1,-1,1}, {-1,-1,1}, {1,-1,-1}, {-1,-1,-1}};

    int * h_sgnv = new int [NSWEEPS * MESHDIM]();
    int * h_sgnt = new int [NSWEEPS * MESHDIM](); 

    for (int index = 0; index < NSWEEPS * MESHDIM; index++)
    {
        int j = index / NSWEEPS;
    	int i = index % NSWEEPS;				

	    h_sgnv[i + j * NSWEEPS] = sgnv[i][j];
	    h_sgnt[i + j * NSWEEPS] = sgnt[i][j];
    }

    T = new float[volsize]();

    cudaMalloc((void**)&(d_T), volsize*sizeof(float));
    cudaMalloc((void**)&(d_S), volsize*sizeof(float));

    cudaMalloc((void**)&(d_sgnv), NSWEEPS*MESHDIM*sizeof(int));
    cudaMalloc((void**)&(d_sgnt), NSWEEPS*MESHDIM*sizeof(int));

    cudaMemcpy(d_S, S, volsize * sizeof(float), cudaMemcpyHostToDevice);    

    cudaMemcpy(d_sgnv, h_sgnv, NSWEEPS*MESHDIM*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sgnt, h_sgnt, NSWEEPS*MESHDIM*sizeof(int), cudaMemcpyHostToDevice);

    std::vector<std::vector<int>>().swap(sgnv);
    std::vector<std::vector<int>>().swap(sgnt);

    delete[] h_sgnt;
    delete[] h_sgnv;
}

void Modeling::initialization()
{
    float sx = geometry->xsrc[geometry->sInd[srcId]]; 
    float sy = geometry->ysrc[geometry->sInd[srcId]]; 
    float sz = geometry->zsrc[geometry->sInd[srcId]]; 

    sIdx = (int)((sx + 0.5f*dx) / dx) + nb;
    sIdy = (int)((sy + 0.5f*dy) / dy) + nb;
    sIdz = (int)((sz + 0.5f*dz) / dz) + nb;

    time_set<<<nBlocks,nThreads>>>(d_T, volsize);

    dim3 grid(1,1,1);
    dim3 block(MESHDIM,MESHDIM,MESHDIM);

    time_init<<<grid,block>>>(d_T,d_S,sx,sy,sz,dx,dy,dz,sIdx,sIdy,sIdz,nxx,nzz,nb);
}

void Modeling::eikonal_solver()
{
    for (int sweep = 0; sweep < NSWEEPS; sweep++)
    { 
	    int start = (sweep == 3 || sweep == 5 || sweep == 6 || sweep == 7) ? total_levels : MESHDIM;
	    int end = (start == MESHDIM) ? total_levels + 1 : MESHDIM - 1;
	    int incr = (start == MESHDIM) ? true : false;

	    int xSweepOff = (sweep == 3 || sweep == 4) ? nxx : 0;
	    int ySweepOff = (sweep == 2 || sweep == 5) ? nyy : 0;
	    int zSweepOff = (sweep == 1 || sweep == 6) ? nzz : 0;
		
	    for (int level = start; level != end; level = (incr) ? level + 1 : level - 1)
	    {			
            int xs = max(1, level - (nyy + nzz));	
            int ys = max(1, level - (nxx + nzz));
            
            int xe = min(nxx, level - (MESHDIM - 1));
            int ye = min(nyy, level - (MESHDIM - 1));	
            
            int xr = xe - xs + 1;
            int yr = ye - ys + 1;

            int nThrds = xr * yr;
                
            dim3 bs(16, 16, 1);

            if (nThrds < 32) { bs.x = xr; bs.y = yr; }  

            dim3 gs(iDivUp(xr, bs.x), iDivUp(yr , bs.y), 1);
                
            int sgni = sweep + 0*NSWEEPS;
            int sgnj = sweep + 1*NSWEEPS;
            int sgnk = sweep + 2*NSWEEPS;

            inner_sweep<<<gs, bs>>>(d_S, d_T, d_sgnt, d_sgnv, sgni, sgnj, sgnk, level, xs, ys, 
                                    xSweepOff, ySweepOff, zSweepOff, nxx, nyy, nzz, dx, dy, dz, 
                                    dx2i, dy2i, dz2i, dz2dx2, dz2dy2, dx2dy2, dsum);
	    }
    }
}

void Modeling::compute_seismogram()
{
    int spread = 0;

    float P[4][4][4];

    cudaMemcpy(T, d_T, volsize*sizeof(float), cudaMemcpyDeviceToHost);    

    for (recId = geometry->iRec[srcId]; recId < geometry->fRec[srcId]; recId++)
    {
        float x = geometry->xrec[recId];
        float y = geometry->yrec[recId];
        float z = geometry->zrec[recId];

        float x0 = floorf((x + 0.5f*dx) / dx) * dx;
        float y0 = floorf((y + 0.5f*dy) / dy) * dy;
        float z0 = floorf((z + 0.5f*dz) / dz) * dz;

        float x1 = floorf((x + 0.5f*dx) / dx) * dx + dx;
        float y1 = floorf((y + 0.5f*dy) / dy) * dy + dy;
        float z1 = floorf((z + 0.5f*dz) / dz) * dz + dz;

        float xd = (x - x0) / (x1 - x0);
        float yd = (y - y0) / (y1 - y0);
        float zd = (z - z0) / (z1 - z0);

        int i = (int)((z + 0.5f*dz) / dz) + nb; 
        int j = (int)((x + 0.5f*dx) / dx) + nb;   
        int k = (int)((y + 0.5f*dy) / dy) + nb;         

        for (int pIdx = 0; pIdx < 4; pIdx++)
        {
            for (int pIdy = 0; pIdy < 4; pIdy++)
            {
                for (int pIdz = 0; pIdz < 4; pIdz++)
                {    
                    P[pIdx][pIdy][pIdz] = T[(i + pIdz - 1) + (j + pIdx - 1)*nzz + (k + pIdy - 1)*nxx*nzz];
                }
            }
        }   

        seismogram[spread++] = cubic3d(P, xd, yd, zd);
    }
}

float Modeling::cubic1d(float P[4], float dx)
{
    return P[1] + 0.5f*dx*(P[2] - P[0] + dx*(2.0f*P[0] - 5.0f*P[1] + 4.0f*P[2] - P[3] + dx*(3.0f*(P[1] - P[2]) + P[3] - P[0])));
}

float Modeling::cubic2d(float P[4][4], float dx, float dy)
{    
    float p[4];
    p[0] = cubic1d(P[0], dy);
    p[1] = cubic1d(P[1], dy);
    p[2] = cubic1d(P[2], dy);
    p[3] = cubic1d(P[3], dy);    
    return cubic1d(p, dx);
}

float Modeling::cubic3d(float P[4][4][4], float dx, float dy, float dz)
{    
    float p[4];
    p[0] = cubic2d(P[0], dy, dz);
    p[1] = cubic2d(P[1], dy, dz);
    p[2] = cubic2d(P[2], dy, dz);
    p[3] = cubic2d(P[3], dy, dz);
    return cubic1d(p, dx);
}

void Modeling::export_seismogram()
{    
    std::string data_file = data_folder + modeling_type + "_nStations" + std::to_string(geometry->spread[srcId]) + "_shot_" + std::to_string(geometry->sInd[srcId]+1) + ".bin";
    export_binary_float(data_file, seismogram, geometry->spread[srcId]);    
}

void Modeling::expand_boundary(float * input, float * output)
{
    # pragma omp parallel for
    for (int index = 0; index < nPoints; index++)
    {
        int k = (int) (index / (nx*nz));         
        int j = (int) (index - k*nx*nz) / nz;    
        int i = (int) (index - j*nz - k*nx*nz);  

        output[(i + nb) + (j + nb)*nzz + (k + nb)*nxx*nzz] = input[i + j*nz + k*nx*nz];       
    }

    for (int k = nb; k < nyy - nb; k++)
    {   
        for (int j = nb; j < nxx - nb; j++)
        {
            for (int i = 0; i < nb; i++)            
            {
                output[i + j*nzz + k*nxx*nzz] = input[0 + (j - nb)*nz + (k - nb)*nx*nz];
                output[(nzz - i - 1) + j*nzz + k*nxx*nzz] = input[(nz - 1) + (j - nb)*nz + (k - nb)*nx*nz];
            }
        }
    }

    for (int k = 0; k < nyy; k++)
    {   
        for (int j = 0; j < nb; j++)
        {
            for (int i = 0; i < nzz; i++)
            {
                output[i + j*nzz + k*nxx*nzz] = output[i + nb*nzz + k*nxx*nzz];
                output[i + (nxx - j - 1)*nzz + k*nxx*nzz] = output[i + (nxx - nb - 1)*nzz + k*nxx*nzz];
            }
        }
    }

    for (int k = 0; k < nb; k++)
    {   
        for (int j = 0; j < nxx; j++)
        {
            for (int i = 0; i < nzz; i++)
            {
                output[i + j*nzz + k*nxx*nzz] = output[i + j*nzz + nb*nxx*nzz];
                output[i + j*nzz + (nyy - k - 1)*nxx*nzz] = output[i + j*nzz + (nyy - nb - 1)*nxx*nzz];
            }
        }
    }
}

void Modeling::reduce_boundary(float * input, float * output)
{
    # pragma omp parallel for
    for (int index = 0; index < nPoints; index++)
    {
        int k = (int) (index / (nx*nz));         
        int j = (int) (index - k*nx*nz) / nz;    
        int i = (int) (index - j*nz - k*nx*nz);  

        output[i + j*nz + k*nx*nz] = input[(i + nb) + (j + nb)*nzz + (k + nb)*nxx*nzz];
    }
}

void Modeling::show_information()
{
    auto clear = system("clear");
    
    std::cout << "-------------------------------------------------------------------------------\n";
    std::cout << "                                 \033[34mSeisFAT3D\033[0;0m\n";
    std::cout << "-------------------------------------------------------------------------------\n\n";

    std::cout << "Model dimensions: (z = " << (nz - 1)*dz << ", x = " << (nx - 1) * dx <<", y = " << (ny - 1) * dy << ") m\n\n";

    std::cout << "Running shot " << srcId + 1 << " of " << geometry->nrel << " in total\n\n";

    std::cout << "Current shot position: (z = " << geometry->zsrc[geometry->sInd[srcId]] << 
                                       ", x = " << geometry->xsrc[geometry->sInd[srcId]] << 
                                       ", y = " << geometry->ysrc[geometry->sInd[srcId]] << ") m\n\n";

    std::cout << modeling_name << "\n";
}

void Modeling::compression(float * input, uintc * output, int volsize, float &max_value, float &min_value)
{
    max_value =-1e20f;
    min_value = 1e20f;
    
    # pragma omp parallel for
    for (int index = 0; index < volsize; index++)
    {
        min_value = std::min(input[index], min_value);
        max_value = std::max(input[index], max_value);        
    }

    # pragma omp parallel for
    for (int index = 0; index < volsize; index++)
        output[index] = static_cast<uintc>(1.0f + (COMPRESS - 1)*(input[index] - min_value) / (max_value - min_value));
}

int Modeling::iDivUp(int a, int b) 
{ 
    return ( (a % b) != 0 ) ? (a / b + 1) : (a / b); 
}

__global__ void time_set(float * T, int volsize)
{
    int index = threadIdx.x + blockIdx.x * blockDim.x;

    if (index < volsize) T[index] = 1e6f;
}

__global__ void time_init(float * T, float * S, float sx, float sy, float sz, float dx, float dy, 
                          float dz, int sIdx, int sIdy, int sIdz, int nxx, int nzz, int nb)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    int k = threadIdx.z + blockIdx.z * blockDim.z;

    int yi = sIdy + (k - 1);
    int xi = sIdx + (j - 1);
    int zi = sIdz + (i - 1);

    int index = zi + xi*nzz + yi*nxx*nzz;

    T[index] = S[index] * sqrtf(powf((xi - nb)*dx - sx, 2.0f) + 
                                powf((yi - nb)*dy - sy, 2.0f) +
                                powf((zi - nb)*dz - sz, 2.0f));
}

__global__ void inner_sweep(float * S, float * T, int * sgnt, int * sgnv, int sgni, int sgnj, int sgnk, 
                            int level, int xOffset, int yOffset, int xSweepOffset, int ySweepOffset, int zSweepOffset, 
                            int nxx, int nyy, int nzz, float dx, float dy, float dz, float dx2i, float dy2i, float dz2i, 
                            float dz2dx2, float dz2dy2, float dx2dy2, float dsum)
{
    int x = (blockIdx.x * blockDim.x + threadIdx.x) + xOffset;
    int y = (blockIdx.y * blockDim.y + threadIdx.y) + yOffset;

    float ta, tb, tc, t1, t2, t3, Sref;
    float t1D1, t1D2, t1D3, t1D, t2D1, t2D2, t2D3, t2D, t3D;

    if ((x < nxx) && (y < nyy)) 
    {
	    int z = level - (x + y);
		
        if ((z >= 0) && (z < nzz))	
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
