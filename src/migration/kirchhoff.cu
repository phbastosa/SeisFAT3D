# include "kirchhoff.cuh"

void Kirchhoff::set_specifications()
{
    cudaMalloc((void**)&(d_image), nPoints*sizeof(float));

    cudaMalloc((void**)&(d_Tr), modeling->nPoints*sizeof(float));
    cudaMalloc((void**)&(d_Ts), modeling->nPoints*sizeof(float));

    cudaMalloc((void**)&(d_seismic), modeling->nt*modeling->max_spread*sizeof(float));

    nThreads = 256;
    nBlocks = (int)((nPoints + nThreads - 1) / nThreads);
}

void Kirchhoff::run_cross_correlation()
{
    cudaMemset(d_image, 0.0f, nPoints*sizeof(float));

    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nrel; modeling->srcId++)
    {
        read_seismic_data();

        modeling->show_information();

        std::cout << "\nKirchhoff depth migration: computing image volume\n";

        modeling->forward_solver();

        modeling->reduce_boundary(modeling->T, Ts);

        cudaMemcpy(d_Ts, Ts, modeling->nPoints*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_seismic, seismic, modeling->nt*modeling->geometry->spread[modeling->srcId]*sizeof(float), cudaMemcpyHostToDevice);

        int spread = 0;

        float sx = modeling->geometry->xsrc[modeling->geometry->sInd[modeling->srcId]];
        float sy = modeling->geometry->ysrc[modeling->geometry->sInd[modeling->srcId]];

        for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
        {
            import_binary_float(output_table_folder + "traveltimes_receiver_" + std::to_string(modeling->recId+1) + ".bin", Tr, modeling->nPoints);

            float rx = modeling->geometry->xrec[modeling->recId];
            float ry = modeling->geometry->yrec[modeling->recId];

            float cmp_x = sx + 0.5f*(rx - sx);
            float cmp_y = sy + 0.5f*(ry - sy);
            
            cudaMemcpy(d_Tr, Tr, modeling->nPoints*sizeof(float), cudaMemcpyHostToDevice);

            cross_correlation<<<nBlocks, nThreads>>>(d_Ts, d_Tr, d_image, d_seismic, aperture_x, aperture_y, cmp_x, cmp_y, spread, nx, ny, nz, dx, dy, dz, 
                                                     modeling->nx, modeling->ny, modeling->nz, modeling->dx, modeling->dy, modeling->dz, scale, modeling->nt, modeling->dt);
        
            ++spread;
        }    
    }

    cudaMemcpy(image, d_image, nPoints*sizeof(float), cudaMemcpyDeviceToHost);

    # pragma omp parallel for
    for (int index = 0; index < nPoints; index++)
        image[index] *= 1.0f / (float)(modeling->geometry->nrel);
}

__global__ void cross_correlation(float * Ts, float * Tr, float * image, float * seismic, float aperture_x, float aperture_y, float cmp_x, float cmp_y, int spread, int nx, int ny, int nz, float dx, float dy, float dz, int snx, int sny, int snz, float sdx, float sdy, float sdz, float scale, int nt, float dt)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    float P[4][4][4];

    int k = (int) (index / (nx*nz));         
    int j = (int) (index - k*nx*nz) / nz;   
    int i = (int) (index - j*nz - k*nx*nz); 

    if ((i > scale) && (i < nz - scale - 1) && (j > scale) && (j < nx - scale - 1) && (k > scale) && (k < ny - scale - 1))
    {
        float z = i*dz;
        float x = j*dx;
        float y = k*dy;

        float x0 = floorf(x / sdx) * sdx;
        float y0 = floorf(y / sdx) * sdy;
        float z0 = floorf(z / sdx) * sdz;

        float x1 = floorf(x / sdx) * sdx + sdx;
        float y1 = floorf(y / sdy) * sdy + sdy;
        float z1 = floorf(z / sdz) * sdz + sdz;

        float xd = (x - x0) / (x1 - x0);
        float yd = (y - y0) / (y1 - y0);
        float zd = (z - z0) / (z1 - z0);        

        int idz = (int)((float)(i) / scale); 
        int idx = (int)((float)(j) / scale); 
        int idy = (int)((float)(k) / scale); 

        for (int pIdx = 0; pIdx < 4; pIdx++)
        {
            for (int pIdy = 0; pIdy < 4; pIdy++)
            {
                for (int pIdz = 0; pIdz < 4; pIdz++)
                {    
                    P[pIdx][pIdy][pIdz] = Ts[(idz + pIdz - 1) + (idx + pIdx - 1)*snz + (idy + pIdy - 1)*snx*snz] + 
                                          Tr[(idz + pIdz - 1) + (idx + pIdx - 1)*snz + (idy + pIdy - 1)*snx*snz];
                }
            }
        }   

        float T = cubic3d(P, xd, yd, zd);

        float T = Ts[index] + Tr[index];

        int tId = (int)(T / dt);

        float sigma_x = tanf(aperture_x * M_PI / 180.0f)*i*dz;        
        float sigma_y = tanf(aperture_y * M_PI / 180.0f)*i*dz;        

        float par_x = powf((j*dx - cmp_x) / (sigma_x + 1e-6f), 2.0f);
        float par_y = powf((k*dy - cmp_y) / (sigma_y + 1e-6f), 2.0f);

        float value = expf(-0.5f*(par_x + par_y));

        if (tId < nt) image[index] += value*seismic[tId + spread*nt];  
    }
}

__device__ float cubic1d(float P[4], float dx)
{
    return P[1] + 0.5f*dx*(P[2] - P[0] + dx*(2.0f*P[0] - 5.0f*P[1] + 4.0f*P[2] - P[3] + dx*(3.0f*(P[1] - P[2]) + P[3] - P[0])));
}

__device__ float cubic2d(float P[4][4], float dx, float dy)
{    
    float p[4];
    p[0] = cubic1d(P[0], dy);
    p[1] = cubic1d(P[1], dy);
    p[2] = cubic1d(P[2], dy);
    p[3] = cubic1d(P[3], dy);    
    return cubic1d(p, dx);
}

__device__ float cubic3d(float P[4][4][4], float dx, float dy, float dz)
{    
    float p[4];
    p[0] = cubic2d(P[0], dy, dz);
    p[1] = cubic2d(P[1], dy, dz);
    p[2] = cubic2d(P[2], dy, dz);
    p[3] = cubic2d(P[3], dy, dz);
    return cubic1d(p, dx);
}
