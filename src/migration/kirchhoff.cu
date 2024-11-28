# include "kirchhoff.cuh"

void Kirchhoff::set_specifications()
{
    cudaMalloc((void**)&(d_Tr), modeling->nPoints*sizeof(float));
    cudaMalloc((void**)&(d_Ts), modeling->nPoints*sizeof(float));

    cudaMalloc((void**)&(d_image), modeling->nPoints*sizeof(float));

    cudaMalloc((void**)&(d_seismic), modeling->nt*modeling->max_spread*sizeof(float));

    nThreads = 256;
    nBlocks = (int)((modeling->nPoints + nThreads - 1) / nThreads);
}

void Kirchhoff::run_cross_correlation()
{
    cudaMemset(d_image, 0.0f, modeling->nPoints*sizeof(float));

    for (modeling->srcId = 0; modeling->srcId < modeling->geometry->nrel; modeling->srcId++)
    {
        read_seismic_data();

        modeling->show_information();

        std::cout << "\nKirchhoff depth migration: computing image volume\n";

        modeling->initialization();
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

            cross_correlation<<<nBlocks, nThreads>>>(d_Ts, d_Tr, d_image, d_seismic, aperture_x, aperture_y, cmp_x, cmp_y, spread, modeling->nx, modeling->ny, modeling->nz, modeling->nt, modeling->dt, modeling->dx, modeling->dy, modeling->dz);
        
            ++spread;
        }
    }

    cudaMemcpy(image, d_image, modeling->nPoints*sizeof(float), cudaMemcpyDeviceToHost);

    # pragma omp parallel for
    for (int index = 0; index < modeling->nPoints; index++)
        image[index] *= 1.0f / modeling->geometry->nrel;
}

__global__ void cross_correlation(float * Ts, float * Tr, float * image, float * seismic, float aperture_x, float aperture_y, float cmp_x, float cmp_y, int spread, int nx, int ny, int nz, int nt, float dt, float dx, float dy, float dz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index < nx*ny*nz)
    {
        int k = (int) (index / (nx*nz));         
        int j = (int) (index - k*nx*nz) / nz;   
        int i = (int) (index - j*nz - k*nx*nz); 

        float sigma_x = tanf(aperture_x * PI / 180.0f)*i*dz;        
        float sigma_y = tanf(aperture_y * PI / 180.0f)*i*dz;        

        float par_x = powf((j*dx - cmp_x) / (sigma_x + 1e-6f), 2.0f);
        float par_y = powf((k*dy - cmp_y) / (sigma_y + 1e-6f), 2.0f);

        float value = expf(-0.5*(par_x + par_y));

        float T = Ts[index] + Tr[index]; 
    
        int tId = (int)(T / dt);

        if (tId < nt) image[index] += value * seismic[tId + spread*nt];  
    }
}
