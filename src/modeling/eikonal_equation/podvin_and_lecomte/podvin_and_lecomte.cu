# include "podvin_and_lecomte.cuh"

void Podvin_and_Lecomte::set_modeling_message()
{
    std::cout<<"Running:\n";
    std::cout<<"[0] - Podvin & Lecomte (1991)\n\n"; 
}

void Podvin_and_Lecomte::set_model_boundaries()
{
    nbxl = 1; nbxr = 1;
    nbyl = 1; nbyr = 1;
    nbzu = 1; nbzd = 1;
}

void Podvin_and_Lecomte::set_preconditioners()
{
    modeling_method = std::string("pod");

    K = new float[volsize]();

    check_spatial_spacing();

    cudaMalloc((void**)&(d_S), volsize*sizeof(float));     
    cudaMalloc((void**)&(d_T), volsize*sizeof(float));     
    cudaMalloc((void**)&(d_K), volsize*sizeof(float));     
    cudaMalloc((void**)&(d_nK), volsize*sizeof(float));     
    cudaMalloc((void**)&(d_nT), volsize*sizeof(float));   

    cudaMemcpy(d_S, S, volsize*sizeof(float), cudaMemcpyHostToDevice);
}

void Podvin_and_Lecomte::check_spatial_spacing()
{
    if ((dx != dy) || (dx != dz))
        throw std::invalid_argument("\033[31mError: For Podvin and Lecomte method the model spacing must to be fixed (dx = dy = dz).\033[0;0m");
}

void Podvin_and_Lecomte::initial_setup()
{
    nit = 0;
    
    int sidx = (int)(geometry->shots.x[shot_id] / dx) + nbxl;
    int sidy = (int)(geometry->shots.y[shot_id] / dy) + nbyl;
    int sidz = (int)(geometry->shots.z[shot_id] / dz) + nbzu;

    int sId = sidz + sidx*nzz + sidy*nxx*nzz;

    float t0 = S[sId] * sqrtf(powf((float)((sidx-nbxl)*dx) - geometry->shots.x[shot_id], 2.0f) +
                              powf((float)((sidy-nbyl)*dy) - geometry->shots.y[shot_id], 2.0f) +
                              powf((float)((sidz-nbzu)*dz) - geometry->shots.z[shot_id], 2.0f));

    for (int index = 0; index < volsize; index++)
    {    
        T[index] = 1e6f;
        K[index] = 0.0f;
    }

    T[sId] = S[sId] * sqrtf(powf((sidx-nbxl)*dx - geometry->shots.x[shot_id], 2.0f) + powf((sidy-nbyl)*dy - geometry->shots.y[shot_id], 2.0f) + powf((sidz-nbzu)*dz - geometry->shots.z[shot_id], 2.0f));

    int aux = (int)sqrtf(powf(sidx, 2.0f) + powf(sidy,2.0f) + powf(sidz,2.0f)); 
    if (aux > nit) nit = aux;

    aux = (int)sqrtf(powf(nxx - sidx,2.0f) + powf(sidy,2.0f) + powf(sidz,2.0f));
    if (aux > nit) nit = aux;

    aux = (int)sqrtf(powf(sidx,2.0f) + powf(nyy - sidy,2.0f) + powf(sidz,2.0f)); 
    if (aux > nit) nit = aux;

    aux = (int)sqrtf(powf(sidx,2.0f) + powf(sidy,2.0f) + powf(nzz - sidz,2.0f)); 
    if (aux > nit) nit = aux;

    aux = (int)sqrtf(powf(sidx,2.0f) + powf(nyy - sidy,2.0f) + powf(nzz - sidz,2.0f));
    if (aux > nit) nit = aux;

    aux = (int)sqrtf(powf(nxx - sidx,2.0f) + powf(sidy,2.0f) + powf(nzz - sidz,2.0f));
    if (aux > nit) nit = aux;

    aux = (int)sqrtf(powf(nxx - sidx,2.0f) + powf(nyy - sidy,2.0f) + powf(sidz,2.0f));
    if (aux > nit) nit = aux;

    aux = (int)sqrtf(powf(nxx - sidx,2.0f) + powf(nyy - sidy,2.0f) + powf(nzz - sidz,2.0f));
    if (aux > nit) nit = aux;

    K[sId - 1] = 1.0f;
    K[sId + 1] = 1.0f;
    K[sId - nzz] = 1.0f;
    K[sId + nzz] = 1.0f;
    K[sId - nxx*nzz] = 1.0f;
    K[sId + nxx*nzz] = 1.0f;
    K[sId + 1 - nzz] = 1.0f;
    K[sId - 1 - nzz] = 1.0f;
    K[sId + 1 + nzz] = 1.0f;
    K[sId - 1 + nzz] = 1.0f;
    K[sId + 1 + nxx*nzz] = 1.0f;
    K[sId + 1 - nxx*nzz] = 1.0f;
    K[sId - 1 + nxx*nzz] = 1.0f;
    K[sId - 1 - nxx*nzz] = 1.0f;
    K[sId - nzz - nxx*nzz] = 1.0f;
    K[sId - nzz + nxx*nzz] = 1.0f;
    K[sId + nzz - nxx*nzz] = 1.0f;
    K[sId + nzz + nxx*nzz] = 1.0f;
    K[sId + 1 + nzz + nxx*nzz] = 1.0f;
    K[sId + 1 + nzz - nxx*nzz] = 1.0f;
    K[sId + 1 - nzz + nxx*nzz] = 1.0f;
    K[sId + 1 - nzz - nxx*nzz] = 1.0f;
    K[sId - 1 - nzz - nxx*nzz] = 1.0f;
    K[sId - 1 - nzz + nxx*nzz] = 1.0f;
    K[sId - 1 + nzz - nxx*nzz] = 1.0f;
    K[sId - 1 + nzz + nxx*nzz] = 1.0f;
}

void Podvin_and_Lecomte::forward_solver()
{
    cudaMemcpy(d_K, K, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_T, T, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_nT, T, volsize*sizeof(float), cudaMemcpyHostToDevice);

    int blocksPerGrid = volsize / threadsPerBlock;

    float dh = dx;

    for (int it = 0; it < nit; it++)
    {
        fdm_operators<<<blocksPerGrid,threadsPerBlock>>>(d_S, d_T, d_K, d_nT, dh, nxx, nyy, nzz);        
        cudaDeviceSynchronize();

        cudaMemset(d_nK, 0.0f, volsize*sizeof(float));

        expanding_box<<<blocksPerGrid,threadsPerBlock>>>(d_K, d_nK, nxx, nyy, nzz);
        cudaDeviceSynchronize();

        update_volume<<<blocksPerGrid,threadsPerBlock>>>(d_T, d_nT, d_K, d_nK, volsize);
        cudaDeviceSynchronize();
    }

    cudaMemcpy(T, d_T, volsize*sizeof(float), cudaMemcpyDeviceToHost);
}

void Podvin_and_Lecomte::free_space()
{
    delete[] K;

    cudaFree(d_K);
    cudaFree(d_S);
    cudaFree(d_T);
    cudaFree(d_nK);
    cudaFree(d_nT);
}

__global__ void fdm_operators(float * S, float * T, float * K, float * nT, float h, int nxx, int nyy, int nzz)
{
    float sqrt2 = sqrtf(2.0f);
    float sqrt3 = sqrtf(3.0f);

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    
    if ((K[index] == 1.0f) && (index < nxx*nyy*nzz))
    {
        int k = (int) (index / (nxx*nzz));         // y direction
        int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
        int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

        if ((i > 0) && (i < nzz-1) && (j > 0) && (j < nxx-1) && (k > 0) && (k < nyy-1))
        {
            float lowest = T[index];
            float Tijk, T1, T2, Sref, M, N, P, Q, hs2;    

            /* 1D operator head wave: i,j-1,k -> i,j,k (x direction) */
            Tijk = T[index - nzz] + h * min(S[index - nzz], 
                                        min(S[index - 1 - nzz], 
                                        min(S[index - nzz - nxx*nzz], S[index - 1 - nzz - nxx*nzz]))); 
            
            if (Tijk < lowest) lowest = Tijk;

            /* 1D operator head wave: i,j+1,k -> i,j,k (x direction) */
            Tijk = T[index + nzz] + h * min(S[index], 
                                        min(S[index - 1], 
                                        min(S[index - nxx*nzz], S[index - 1 - nxx*nzz])));
            
            if (Tijk < lowest) lowest = Tijk;

            /* 1D operator head wave: i,j,k-1 -> i,j,k (y direction) */
            Tijk = T[index - nxx*nzz] + h * min(S[index - nxx*nzz], 
                                            min(S[index - nzz - nxx*nzz], 
                                            min(S[index - 1 - nxx*nzz], S[index - 1 - nzz - nxx*nzz]))); 
            
            if (Tijk < lowest) lowest = Tijk;

            /* 1D operator head wave: i,j,k+1 -> i,j,k (y direction) */
            Tijk = T[index + nxx*nzz] + h * min(S[index],
                                            min(S[index - 1], 
                                            min(S[index - nzz], S[index - 1 - nzz]))); 
            
            if (Tijk < lowest) lowest = Tijk;

            /* 1D operator head wave: i-1,j,k -> i,j,k (z direction) */
            Tijk = T[index - 1] + h * min(S[index - 1], 
                                      min(S[index - 1 - nzz], 
                                      min(S[index - 1 - nxx*nzz], S[index - 1 - nzz - nxx*nzz]))); 
            
            if (Tijk < lowest) lowest = Tijk;

            /* 1D operator head wave: i+1,j,k -> i,j,k (z direction) */
            Tijk = T[index + 1] + h * min(S[index], 
                                      min(S[index - nzz], 
                                      min(S[index - nxx*nzz], S[index - nzz - nxx*nzz]))); 
            
            if (Tijk < lowest) lowest = Tijk;
                
            /* 1D operator diffraction XZ plane */
            
            // i-1,j-1,k -> i,j,k
            Tijk = T[index - 1 - nzz] + h*sqrt2*min(S[index - 1 - nzz], S[index - 1 - nzz - nxx*nzz]); 
            if (Tijk < lowest) lowest = Tijk;

            // i-1,j+1,k -> i,j,k
            Tijk = T[index - 1 + nzz] + h*sqrt2*min(S[index - 1], S[index - 1 - nxx*nzz]); 
            if (Tijk < lowest) lowest = Tijk;
            
            // i+1,j-1,k -> i,j,k
            Tijk = T[index + 1 - nzz] + h*sqrt2*min(S[index - nzz], S[index - nzz - nxx*nzz]); 
            if (Tijk < lowest) lowest = Tijk;
            
            // i+1,j+1,k -> i,j,k
            Tijk = T[index + 1 + nzz] + h*sqrt2*min(S[index], S[index - nxx*nzz]); 
            if (Tijk < lowest) lowest = Tijk;

            /* 1D operator diffraction YZ plane */

            // i-1,j,k-1 -> i,j,k
            Tijk = T[index - 1 - nxx*nzz] + h*sqrt2*min(S[index - 1 - nxx*nzz], S[index - 1 - nzz - nxx*nzz]); 
            if (Tijk < lowest) lowest = Tijk;

            // i-1,j,k+1 -> i,j,k
            Tijk = T[index - 1 + nxx*nzz] + h*sqrt2*min(S[index - 1], S[index - 1 - nzz]); 
            if (Tijk < lowest) lowest = Tijk;
            
            // i+1,j,k-1 -> i,j,k
            Tijk = T[index + 1 - nxx*nzz] + h*sqrt2*min(S[index - nxx*nzz], S[index - nzz - nxx*nzz]); 
            if (Tijk < lowest) lowest = Tijk;
            
            // i+1,j,k+1 -> i,j,k
            Tijk = T[index + 1 + nxx*nzz] + h*sqrt2*min(S[index], S[index - nzz]); 
            if (Tijk < lowest) lowest = Tijk;

            /* 1D operator diffraction XY plane */
            
            // i,j-1,k-1 -> i,j,k
            Tijk = T[index - nzz - nxx*nzz] + h*sqrt2*min(S[index - nzz - nxx*nzz], S[index - 1 - nzz - nxx*nzz]); 
            if (Tijk < lowest) lowest = Tijk;

            // i,j-1,k+1 -> i,j,k
            Tijk = T[index - nzz + nxx*nzz] + h*sqrt2*min(S[index - nzz], S[index - 1 - nzz]); 
            if (Tijk < lowest) lowest = Tijk;

            // i,j+1,k-1 -> i,j,k
            Tijk = T[index + nzz - nxx*nzz] + h*sqrt2*min(S[index - nxx*nzz], S[index - 1 - nxx*nzz]); 
            if (Tijk < lowest) lowest = Tijk;

            // i,j+1,k+1 -> i,j,k
            Tijk = T[index + nzz + nxx*nzz] + h*sqrt2*min(S[index], S[index - 1]); 
            if (Tijk < lowest) lowest = Tijk;

            /* 1D operator corner diffractions */

            // i-1,j-1,k-1 -> i,j,k
            Tijk = T[index - 1 - nzz - nxx*nzz] + h*sqrt3*S[index - 1 - nzz - nxx*nzz]; 
            if (Tijk < lowest) lowest = Tijk;

            // i-1,j-1,k+1 -> i,j,k
            Tijk = T[index - 1 - nzz + nxx*nzz] + h*sqrt3*S[index - 1 - nzz]; 
            if (Tijk < lowest) lowest = Tijk;

            // i+1,j-1,k-1 -> i,j,k
            Tijk = T[index + 1 - nzz - nxx*nzz] + h*sqrt3*S[index - nzz - nxx*nzz]; 
            if (Tijk < lowest) lowest = Tijk;

            // i+1,j-1,k+1 -> i,j,k
            Tijk = T[index + 1 - nzz + nxx*nzz] + h*sqrt3*S[index - nzz]; 
            if (Tijk < lowest) lowest = Tijk;

            // i-1,j+1,k-1 -> i,j,k
            Tijk = T[index - 1 + nzz - nxx*nzz] + h*sqrt3*S[index - 1 - nxx*nzz]; 
            if (Tijk < lowest) lowest = Tijk;

            // i-1,j+1,k+1 -> i,j,k
            Tijk = T[index - 1 + nzz + nxx*nzz] + h*sqrt3*S[index - 1]; 
            if (Tijk < lowest) lowest = Tijk;

            // i+1,j+1,k-1 -> i,j,k
            Tijk = T[index + 1 + nzz - nxx*nzz] + h*sqrt3*S[index - nxx*nzz]; 
            if (Tijk < lowest) lowest = Tijk;

            // i+1,j+1,k+1 -> i,j,k
            Tijk = T[index + 1 + nzz + nxx*nzz] + h*sqrt3*S[index]; 
            if (Tijk < lowest) lowest = Tijk;

            /* 2D operator XZ plane: First Quadrant*/

            Sref = min(S[index - 1 - nzz], S[index - 1 - nzz - nxx*nzz]);

            // i,j-1,k - i-1,j-1,k -> i,j,k
            T1 = T[index - nzz];
            T2 = T[index - 1 - nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            // i-1,j,k - i-1,j-1,k -> i,j,k
            T1 = T[index - 1];
            T2 = T[index - 1 - nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            /* 2D operator XZ plane: Second Quadrant*/                        

            Sref = min(S[index - nzz], S[index - nzz - nxx*nzz]);

            // i,j-1,k - i+1,j-1,k -> i,j,k
            T1 = T[index - nzz];
            T2 = T[index + 1 - nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            // i+1,j,k - i+1,j-1,k -> i,j,k
            T1 = T[index + 1];
            T2 = T[index + 1 - nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            /* 2D operator XZ plane: Third Quadrant*/                        

            Sref = min(S[index], S[index - nxx*nzz]);

            // i+1,j,k - i+1,j+1,k -> i,j,k
            T1 = T[index + 1];
            T2 = T[index + 1 + nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            // i,j+1,k - i+1,j+1,k -> i,j,k
            T1 = T[index + nzz];
            T2 = T[index + 1 + nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            /* 2D operator XZ plane: Fourth Quadrant*/                        

            Sref = min(S[index - 1], S[index - 1 - nxx*nzz]);

            // i,j+1,k - i-1,j+1,k -> i,j,k
            T1 = T[index + nzz];
            T2 = T[index - 1 + nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            // i-1,j,k - i-1,j+1,k -> i,j,k
            T1 = T[index - 1];
            T2 = T[index - 1 + nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            /* 2D operator YZ plane: First Quadrant */                        

            Sref = min(S[index - 1 - nxx*nzz], S[index - 1 - nzz - nxx*nzz]);

            // i,j,k-1 - i-1,j,k-1 -> i,j,k
            T1 = T[index - nxx*nzz];
            T2 = T[index - 1 - nxx*nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            // i-1,j,k - i-1,j,k-1 -> i,j,k
            T1 = T[index - 1];
            T2 = T[index - 1 - nxx*nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            /* 2D operator YZ plane: Second Quadrant */                        

            Sref = min(S[index - nxx*nzz], S[index - nzz - nxx*nzz]);

            // i,j,k-1 - i+1,j,k-1 -> i,j,k
            T1 = T[index - nxx*nzz];
            T2 = T[index + 1 - nxx*nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            // i+1,j,k - i+1,j,k-1 -> i,j,k
            T1 = T[index + 1];
            T2 = T[index + 1 - nxx*nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            /* 2D operator YZ plane: Third Quadrant*/                        

            Sref = min(S[index], S[index - nzz]);

            // i+1,j,k - i+1,j,k+1 -> i,j,k
            T1 = T[index + 1];
            T2 = T[index + 1 + nxx*nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            // i,j,k+1 - i+1,j,k+1 -> i,j,k
            T1 = T[index + nxx*nzz];
            T2 = T[index + 1 + nxx*nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            /* 2D operator YZ plane: Fourth Quadrant*/                        

            Sref = min(S[index - 1], S[index - 1 - nzz]);

            // i,j,k+1 - i-1,j,k+1 -> i,j,k
            T1 = T[index + nxx*nzz];
            T2 = T[index - 1 + nxx*nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            // i-1,j,k - i-1,j,k+1 -> i,j,k
            T1 = T[index - 1];
            T2 = T[index - 1 + nxx*nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            /* 2D operator XY plane: First Quadrant*/                        

            Sref = min(S[index - nzz - nxx*nzz], S[index - 1 - nzz - nxx*nzz]);

            // i,j-1,k - i,j-1,k-1 -> i,j,k
            T1 = T[index - nzz];
            T2 = T[index - nzz - nxx*nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            // i,j,k-1 - i,j-1,k-1 -> i,j,k
            T1 = T[index - nxx*nzz];
            T2 = T[index - nzz - nxx*nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            /* 2D operator XY plane: Second Quadrant*/                        

            Sref = min(S[index - nzz], S[index - 1 - nzz]);

            // i,j-1,k - i,j-1,k+1 -> i,j,k
            T1 = T[index - nzz];
            T2 = T[index - nzz + nxx*nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            // i,j,k+1 - i,j-1,k+1 -> i,j,k
            T1 = T[index + nxx*nzz];
            T2 = T[index - nzz + nxx*nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            /* 2D operator XY plane: Third Quadrant*/                        

            Sref = min(S[index], S[index - 1]);

            // i,j,k+1 - i,j+1,k+1 -> i,j,k
            T1 = T[index + nxx*nzz];
            T2 = T[index + nzz + nxx*nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            // i,j+1,k - i,j+1,k+1 -> i,j,k
            T1 = T[index + nzz];
            T2 = T[index + nzz + nxx*nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            /* 2D operator XY plane: Fourth Quadrant*/                        

            Sref = min(S[index - nxx*nzz], S[index - 1 - nxx*nzz]);

            // i,j+1,k - i,j+1,k-1 -> i,j,k
            T1 = T[index + nzz];
            T2 = T[index + nzz - nxx*nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            // i,j,k-1 - i,j+1,k-1 -> i,j,k
            T1 = T[index - nxx*nzz];
            T2 = T[index + nzz - nxx*nzz];
            if ((T1 - T2) > 0.0f)
            {
                if ((T1 - T2) < h*Sref/sqrt2)
                {
                    Tijk = T1 + sqrtf(h*h*Sref*Sref - (T1 - T2)*(T1 - T2));
                    if (Tijk < lowest) lowest = Tijk;
                }
            }

            /* 3D operator - First octant: XY plane */

            Sref = S[index - 1 - nzz - nxx*nzz];
            hs2 = h*h*Sref*Sref;

            M = T[index - 1 - nzz - nxx*nzz]; /* i-1,j-1,k-1 */  
            N = T[index - 1 - nzz];           /* i-1,j-1, k  */  
            P = T[index - 1 - nxx*nzz];       /* i-1, j ,k-1 */
            Q = T[index - 1];                 /* i-1, j , k  */

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - First octant: YZ plane */

            M = T[index - 1 - nzz - nxx*nzz]; /* i-1,j-1,k-1 */   
            N = T[index - 1 - nzz];           /* i-1,j-1, k  */   
            P = T[index - nzz - nxx*nzz];     /*  i ,j-1,k-1 */   
            Q = T[index - nzz];               /*  i ,j-1, k  */   

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - First octant: XZ plane */

            M = T[index - 1 - nzz - nxx*nzz]; /* i-1,j-1,k-1 */   
            N = T[index - nzz - nxx*nzz];     /*  i ,j-1,k-1 */         
            P = T[index - 1 - nxx*nzz];       /* i-1, j ,k-1 */ 
            Q = T[index - nxx*nzz];           /*  i , j ,k-1 */       

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
            ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
            ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
            ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
            ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Second octant: XY plane */

            Sref = S[index - 1 - nxx*nzz];
            hs2 = h*h*Sref*Sref;

            M = T[index - 1 + nzz - nxx*nzz];  /* i-1,j+1,k-1 */ 
            N = T[index - 1 - nxx*nzz];        /* i-1, j ,k-1 */     
            P = T[index - 1 + nzz];            /* i-1,j+1, k  */
            Q = T[index - 1];                  /* i-1, j , k  */

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Second octant: YZ plane */

            M = T[index - 1 + nzz - nxx*nzz]; /* i-1,j+1,k-1 */   
            N = T[index - 1 + nzz];           /* i-1,j+1, k  */   
            P = T[index + nzz - nxx*nzz];     /*  i ,j+1,k-1 */   
            Q = T[index + nzz];               /*  i ,j+1, k  */   

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Second octant: XZ plane */

            M = T[index - 1 + nzz - nxx*nzz]; /* i-1,j+1,k-1 */   
            N = T[index - 1 - nxx*nzz];       /* i-1, j ,k-1 */       
            P = T[index + nzz - nxx*nzz];     /*  i ,j+1,k-1 */   
            Q = T[index - nxx*nzz];           /*  i , j ,k-1 */       

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Third octant: XY plane */

            Sref = S[index - 1];
            hs2 = h*h*Sref*Sref;

            M = T[index - 1 + nzz + nxx*nzz]; /* i-1,j+1,k+1 */   
            N = T[index - 1 + nzz];           /* i-1,j+1, k  */   
            P = T[index - 1 + nxx*nzz];       /* i-1, j ,k+1 */ 
            Q = T[index - 1];                 /* i-1, j , k  */ 

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Third octant: YZ plane */

            M = T[index - 1 + nzz + nxx*nzz]; /* i-1,j+1,k+1 */   
            N = T[index + nzz + nxx*nzz];     /*  i ,j+1,k+1 */         
            P = T[index - 1 + nzz];           /* i-1,j+1, k  */ 
            Q = T[index + nzz];               /*  i ,j+1, k  */   

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Third octant: XZ plane */

            M = T[index - 1 + nzz + nxx*nzz]; /* i-1,j+1,k+1 */   
            N = T[index - 1 + nxx*nzz];       /* i-1, j ,k+1 */       
            P = T[index + nzz + nxx*nzz];     /*  i ,j+1,k+1 */   
            Q = T[index + nxx*nzz];           /*  i , j ,k+1 */       

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Fourth octant: XY plane */

            Sref = S[index - 1 - nzz];
            hs2 = h*h*Sref*Sref;

            M = T[index - 1 - nzz + nxx*nzz]; /* i-1,j-1,k+1 */  
            N = T[index - 1 + nxx*nzz];       /* i-1, j ,k+1 */      
            P = T[index - 1 - nzz];           /* i-1,j-1, k  */
            Q = T[index - 1];                 /* i-1, j , k  */

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Fourth octant: YZ plane */

            M = T[index - 1 - nzz + nxx*nzz]; /* i-1,j-1,k+1 */  
            N = T[index - 1 - nzz];           /* i-1,j-1, k  */  
            P = T[index - nzz + nxx*nzz];     /*  i ,j-1,k+1 */  
            Q = T[index - nzz];               /*  i ,j-1, k  */  

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Fourth octant: XZ plane */

            M = T[index - 1 - nzz + nxx*nzz]; /* i-1,j-1,k+1 */   
            N = T[index - nzz + nxx*nzz];     /*  i ,j-1,k+1 */         
            P = T[index - 1 + nxx*nzz];       /* i-1, j ,k+1 */ 
            Q = T[index + nxx*nzz];           /*  i , j ,k+1 */       

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Fifth octant: XY plane */

            Sref = S[index - nzz - nxx*nzz];
            hs2 = h*h*Sref*Sref;

            M = T[index + 1 - nzz - nxx*nzz]; /* i+1,j-1,k-1 */  
            N = T[index + 1 - nxx*nzz];       /* i+1, j ,k-1 */      
            P = T[index + 1 - nzz];           /* i+1,j-1, k  */
            Q = T[index + 1];                 /* i+1, j , k  */

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Fifth octant: YZ plane */

            M = T[index + 1 - nzz - nxx*nzz]; /* i+1,j-1,k-1 */   
            N = T[index + 1 - nzz];           /* i+1,j-1, k  */   
            P = T[index - nzz - nxx*nzz];     /*  i ,j-1,k-1 */   
            Q = T[index - nzz];               /*  i ,j-1, k  */   

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
            ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
            ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
            ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
            ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Fifth octant: XZ plane */

            M = T[index + 1 - nzz - nxx*nzz]; /* i+1,j-1,k-1 */   
            N = T[index - nzz - nxx*nzz];     /*  i ,j-1,k-1 */         
            P = T[index + 1 - nxx*nzz];       /* i+1, j ,k-1 */ 
            Q = T[index - nxx*nzz];           /*  i , j ,k-1 */       

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Sixth octant: XY plane */

            Sref = S[index - nxx*nzz];
            hs2 = h*h*Sref*Sref;

            M = T[index + 1 + nzz - nxx*nzz]; /* i+1,j+1,k-1 */   
            N = T[index + 1 + nzz];           /* i+1,j+1, k  */   
            P = T[index + 1 - nxx*nzz];       /* i+1, j ,k-1 */ 
            Q = T[index + 1];                 /* i+1, j , k  */ 

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Sixth octant: YZ plane */

            M = T[index + 1 + nzz - nxx*nzz]; /* i+1,j+1,k-1 */   
            N = T[index + nzz - nxx*nzz];     /*  i ,j+1,k-1 */        
            P = T[index + 1 + nzz];           /* i+1,j+1, k  */
            Q = T[index + nzz];               /*  i ,j+1, k  */  

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Sixth octant: XZ plane */

            M = T[index + 1 + nzz - nxx*nzz]; /* i+1,j+1,k-1 */   
            N = T[index + 1 - nxx*nzz];       /* i+1, j ,k-1 */       
            P = T[index + nzz - nxx*nzz];     /*  i ,j+1,k-1 */   
            Q = T[index - nxx*nzz];           /*  i , j ,k-1 */       

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Seventh octant: XY plane */
            
            Sref = S[index - nzz];
            hs2 = h*h*Sref*Sref;

            M = T[index + 1 - nzz + nxx*nzz]; /* i+1,j-1,k+1 */   
            N = T[index + 1 - nzz];           /* i+1,j-1, k  */   
            P = T[index + 1 + nxx*nzz];       /* i+1, j ,k+1 */ 
            Q = T[index + 1];                 /* i+1, j , k  */ 

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Seventh octant: YZ plane */

            M = T[index + 1 - nzz + nxx*nzz]; /* i+1,j-1,k+1 */  
            N = T[index - nzz + nxx*nzz];     /*  i ,j-1,k+1 */        
            P = T[index + 1 - nzz];           /* i+1,j-1, k  */
            Q = T[index - nzz];               /*  i ,j-1, k  */  

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Seventh octant: XZ plane */

            M = T[index + 1 - nzz + nxx*nzz]; /* i+1,j-1,k+1 */    
            N = T[index + 1 + nxx*nzz];       /* i+1, j ,k+1 */        
            P = T[index - nzz + nxx*nzz];     /*  i ,j-1,k+1 */    
            Q = T[index + nxx*nzz];           /*  i , j ,k+1 */        

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Eighth octant: XY plane */

            Sref = S[index];
            hs2 = h*h*Sref*Sref;

            M = T[index + 1 + nzz + nxx*nzz]; /* i+1,j+1,k+1 */  
            N = T[index + 1 + nxx*nzz];       /* i+1, j ,k+1 */      
            P = T[index + 1 + nzz];           /* i+1,j+1, k  */
            Q = T[index + 1];                 /* i+1, j , k  */

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Eighth octant: YZ plane */

            M = T[index + 1 + nzz + nxx*nzz]; /* i+1,j+1,k+1 */   
            N = T[index + 1 + nzz];           /* i+1,j+1, k  */   
            P = T[index + nzz + nxx*nzz];     /*  i ,j+1,k+1 */   
            Q = T[index + nzz];               /*  i ,j+1, k  */   

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* 3D operator - Eighth octant: XZ plane */

            M = T[index + 1 + nzz + nxx*nzz]; /* i+1,j+1,k+1 */   
            N = T[index + nzz + nxx*nzz];     /*  i ,j+1,k+1 */         
            P = T[index + 1 + nxx*nzz];       /* i+1, j ,k+1 */ 
            Q = T[index + nxx*nzz];           /*  i , j ,k+1 */       

            // MNP -> R 
            if ((M <= N) && (M <= P) && 
                ((2.0f*(P-M)*(P-M) + (N-M)*(N-M)) <= hs2) && 
                ((2.0f*(N-M)*(N-M) + (P-M)*(P-M)) <= hs2) && 
                ((N-M)*(N-M) + (P-M)*(P-M) + (N-M)*(P-M) >= 0.5f*hs2))
            {
                Tijk = N + P - M + sqrtf(hs2 - (N-M)*(N-M) - (P-M)*(P-M));
                if (Tijk < lowest) lowest = Tijk;
            }   

            // QNP -> R    
            if ((N <= Q) && (P <= Q) && 
                ((Q-N)*(Q-N) + (Q-P)*(Q-P) + (Q-N)*(Q-P) <= 0.5f*hs2))    
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (Q-P)*(Q-P));    
                if (Tijk < lowest) lowest = Tijk;
            }

            // NMQ -> R
            if ((N-M >= 0) && (N-M <= Q-N) && 
                (2*(Q-N)*(Q-N) + (N-M)*(N-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-N)*(Q-N) - (N-M)*(N-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            // PMQ -> R
            if ((P-M >= 0) && (P-M <= Q-P) && 
                (2*(Q-P)*(Q-P) + (P-M)*(P-M) <= hs2))
            {
                Tijk = Q + sqrtf(hs2 - (Q-P)*(Q-P) - (P-M)*(P-M));    
                if (Tijk < lowest) lowest = Tijk;
            }        

            /* Time atualization */
            if (lowest == T[index]) K[index] = 0.0f;

            nT[index] = lowest;
        }
    }
}

__global__ void expanding_box(float * K, float * nK, int nxx, int nyy, int nzz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if ((K[index] == 1.0f) && (index < nxx*nyy*nzz))
    {
        int k = (int) (index / (nxx*nzz));         // y direction
        int j = (int) (index - k*nxx*nzz) / nzz;   // x direction
        int i = (int) (index - j*nzz - k*nxx*nzz); // z direction

        if ((i > 0) && (i < nzz-1) && (j > 0) && (j < nxx-1) && (k > 0) && (k < nyy-1))
        {
            nK[index - 1] = 1.0f;
            nK[index + 1] = 1.0f;
            nK[index - nzz] = 1.0f;
            nK[index + nzz] = 1.0f;
            nK[index - nxx*nzz] = 1.0f;
            nK[index + nxx*nzz] = 1.0f;
            nK[index + 1 - nzz] = 1.0f;
            nK[index - 1 - nzz] = 1.0f;
            nK[index + 1 + nzz] = 1.0f;
            nK[index - 1 + nzz] = 1.0f;
            nK[index + 1 + nxx*nzz] = 1.0f;
            nK[index + 1 - nxx*nzz] = 1.0f;
            nK[index - 1 + nxx*nzz] = 1.0f;
            nK[index - 1 - nxx*nzz] = 1.0f;
            nK[index - nzz - nxx*nzz] = 1.0f;
            nK[index - nzz + nxx*nzz] = 1.0f;
            nK[index + nzz - nxx*nzz] = 1.0f;
            nK[index + nzz + nxx*nzz] = 1.0f;
            nK[index + 1 + nzz + nxx*nzz] = 1.0f;
            nK[index + 1 + nzz - nxx*nzz] = 1.0f;
            nK[index + 1 - nzz + nxx*nzz] = 1.0f;
            nK[index + 1 - nzz - nxx*nzz] = 1.0f;
            nK[index - 1 - nzz - nxx*nzz] = 1.0f;
            nK[index - 1 - nzz + nxx*nzz] = 1.0f;
            nK[index - 1 + nzz - nxx*nzz] = 1.0f;
            nK[index - 1 + nzz + nxx*nzz] = 1.0f;
        }
    }
}

__global__ void update_volume(float * T, float * nT, float * K, float * nK, int N)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index < N)
    {
        T[index] = nT[index];
        K[index] = nK[index];
    }
}