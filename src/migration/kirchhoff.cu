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

        std::cout << "\nKirchhoff depth migration: computing image volume\n\n";

        modeling->initialization();
        modeling->forward_solver();

        modeling->reduce_boundary(modeling->T, Ts);

        cudaMemcpy(d_Ts, Ts, modeling->nPoints*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_seismic, seismic, modeling->nt*modeling->geometry->spread[modeling->srcId]*sizeof(float), cudaMemcpyHostToDevice);

        int spread = 0;

        for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
        {
            import_binary_float("../outputs/travelTimeTables/receiver_travelTimes_" + std::to_string(modeling->recId+1) + ".bin", Tr, modeling->nPoints);
            
            cudaMemcpy(d_Tr, Tr, modeling->nPoints*sizeof(float), cudaMemcpyHostToDevice);

            cross_correlation<<<nBlocks, nThreads>>>(d_seismic, d_Ts, d_Tr, d_image, modeling->nPoints, spread, modeling->nt, modeling->dt);
        
            ++spread;
        }
    }

    cudaMemcpy(image, d_image, modeling->nPoints*sizeof(float), cudaMemcpyDeviceToHost);
}

__global__ void cross_correlation(float * seismic, float * Ts, float * Tr, float * image, int nPoints, int spread, int nt, float dt)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index < nPoints)
    {
        float T = Ts[index] + Tr[index];

        int tId = (int)(T / dt);

        if (tId < nt) image[index] += seismic[tId + spread*nt];
    }
}
