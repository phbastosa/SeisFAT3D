# include "elastic_vti.cuh"

void Elastic_VTI::set_conditions()
{
    modeling_type = "elastic_vti";
    modeling_name = "Modeling type: Elastic vertically transverse isotropic time propagation";

    eikonal = new Eikonal_VTI();
    eikonal->parameters = parameters;
    eikonal->set_parameters();

    T = new float[nPoints]();
    TT = new float[volsize]();

    std::string e_file = catch_parameter("epsilon_model_file", parameters);
    std::string d_file = catch_parameter("delta_model_file", parameters);
    std::string g_file = catch_parameter("gamma_model_file", parameters);

    float * e = new float[nPoints]();
    float * d = new float[nPoints]();
    float * g = new float[nPoints]();

    E = new float[volsize]();
    D = new float[volsize]();
    G = new float[volsize]();

    import_binary_float(e_file, e, nPoints);
    import_binary_float(d_file, d, nPoints);
    import_binary_float(g_file, g, nPoints);

    expand_boundary(e, E);
    expand_boundary(d, D);
    expand_boundary(g, G);

    delete[] e;
    delete[] d;
    delete[] g;

    B = new float[volsize]();
    C11 = new float[volsize]();
    C33 = new float[volsize]();
    C55 = new float[volsize]();
    C66 = new float[volsize]();
    C13 = new float[volsize]();

    # pragma omp parallel for
    for (int index = 0; index < volsize; index++)
    {
        B[index] = 1.0f / Ro[index];
        
        C33[index] = Vp[index]*Vp[index]*Ro[index];
        C55[index] = Vs[index]*Vs[index]*Ro[index];

        C11[index] = C33[index]*(1.0f + 2.0f*E[index]); 
        C66[index] = C55[index]*(1.0f + 2.0f*G[index]);

        C13[index] = sqrtf(2.0f*D[index]*C33[index]*(C33[index] - C55[index]) + powf(C33[index] - C55[index], 2.0f)) - C55[index];
    }

    cudaMalloc((void**)&(d_T), volsize*sizeof(float));

    cudaMalloc((void**)&(d_B), volsize*sizeof(float));
    cudaMalloc((void**)&(d_C11), volsize*sizeof(float));
    cudaMalloc((void**)&(d_C33), volsize*sizeof(float));
    cudaMalloc((void**)&(d_C55), volsize*sizeof(float));
    cudaMalloc((void**)&(d_C66), volsize*sizeof(float));
    cudaMalloc((void**)&(d_C13), volsize*sizeof(float));

    cudaMemcpy(d_B, B, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_C11, C11, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_C33, C33, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_C55, C55, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_C66, C66, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_C13, C13, volsize*sizeof(float), cudaMemcpyHostToDevice);
}

void Elastic_VTI::forward_solver()
{
    eikonal->srcId = srcId;

    eikonal->initialization();
    eikonal->forward_solver();

    eikonal->reduce_boundary(eikonal->T, T);
    
    expand_boundary(T, TT);
    
    cudaMemcpy(d_T, TT, volsize*sizeof(float), cudaMemcpyHostToDevice);

    for (int tId = 0; tId < nt + tlag; tId++)
    {
        compute_pressure<<<nBlocks, nThreads>>>(d_Vx, d_Vy, d_Vz, d_Txx, d_Tyy, d_Tzz, d_Txz, d_Tyz, d_Txy, d_P, d_T, d_C11, d_C33, d_C55, d_C66, d_C13, wavelet, sIdx, sIdy, sIdz, tId, tlag, nt, dx, dy, dz, dt, nxx, nyy, nzz);
        cudaDeviceSynchronize();

        compute_velocity<<<nBlocks, nThreads>>>(d_Vx, d_Vy, d_Vz, d_Txx, d_Tyy, d_Tzz, d_Txz, d_Tyz, d_Txy, d_B, d_T, d1D, d2D, d3D, dx, dy, dz, dt, tId, tlag, nxx, nyy, nzz, nb);
        cudaDeviceSynchronize();

        compute_seismogram<<<sBlocks, nThreads>>>(d_P, rIdx, rIdy, rIdz, seismogram, geometry->spread[srcId], tId, tlag, nt, nxx, nzz);     
        cudaDeviceSynchronize();
    }

    cudaMemcpy(synthetic_data, seismogram, nt*geometry->spread[srcId]*sizeof(float), cudaMemcpyDeviceToHost);
}

__global__ void compute_pressure(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * P, float * T, 
                                 float * C11, float * C33, float * C55, float * C66, float * C13, float * wavelet, int sIdx, int sIdy, int sIdz, 
                                 int tId, int tlag, int nt, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         
    int j = (int) (index - k*nxx*nzz) / nzz;   
    int i = (int) (index - j*nzz - k*nxx*nzz); 

    if ((index == 0) && (tId < nt))
    {
        Txx[sIdz + sIdx*nzz + sIdy*nxx*nzz] += wavelet[tId] / (dx*dy*dz);
        Tyy[sIdz + sIdx*nzz + sIdy*nxx*nzz] += wavelet[tId] / (dx*dy*dz);
        Tzz[sIdz + sIdx*nzz + sIdy*nxx*nzz] += wavelet[tId] / (dx*dy*dz);
    }

    if ((index < nxx*nyy*nzz) && (T[index] < (float)(tId + tlag)*dt))
    {
        if((i >= 3) && (i < nzz-4) && (j >= 3) && (j < nxx-4) && (k >= 3) && (k < nyy-4)) 
        {            
            float dVx_dx = (75.0f*(Vx[i + (j-3)*nzz + k*nxx*nzz] - Vx[i + (j+4)*nzz + k*nxx*nzz]) +
                        1029.0f*(Vx[i + (j+3)*nzz + k*nxx*nzz] - Vx[i + (j-2)*nzz + k*nxx*nzz]) +
                        8575.0f*(Vx[i + (j-1)*nzz + k*nxx*nzz] - Vx[i + (j+2)*nzz + k*nxx*nzz]) +
                        128625.0f*(Vx[i + (j+1)*nzz + k*nxx*nzz] - Vx[i + j*nzz + k*nxx*nzz])) / (107520.0f*dx);

            float dVy_dy = (75.0f*(Vy[i + j*nzz + (k-3)*nxx*nzz] - Vy[i + j*nzz + (k+4)*nxx*nzz]) +
                        1029.0f*(Vy[i + j*nzz + (k+3)*nxx*nzz] - Vy[i + j*nzz + (k-2)*nxx*nzz]) +
                        8575.0f*(Vy[i + j*nzz + (k-1)*nxx*nzz] - Vy[i + j*nzz + (k+2)*nxx*nzz]) +
                        128625.0f*(Vy[i + j*nzz + (k+1)*nxx*nzz] - Vy[i + j*nzz + k*nxx*nzz])) / (107520.0f*dy);

            float dVz_dz = (75.0f*(Vz[(i-3) + j*nzz + k*nxx*nzz] - Vz[(i+4) + j*nzz + k*nxx*nzz]) +
                        1029.0f*(Vz[(i+3) + j*nzz + k*nxx*nzz] - Vz[(i-2) + j*nzz + k*nxx*nzz]) +
                        8575.0f*(Vz[(i-1) + j*nzz + k*nxx*nzz] - Vz[(i+2) + j*nzz + k*nxx*nzz]) +
                        128625.0f*(Vz[(i+1) + j*nzz + k*nxx*nzz] - Vz[i + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

            Txx[index] += dt*(C11[index]*dVx_dx + (C11[index] - 2.0f*C66[index])*dVy_dy + C13[index]*dVz_dz);
            Tyy[index] += dt*(C11[index]*dVy_dy + (C11[index] - 2.0f*C66[index])*dVx_dx + C13[index]*dVz_dz);
            Tzz[index] += dt*(C33[index]*dVz_dz + C13[index]*(dVx_dx + dVy_dy));                    
        }

        if((i >= 0) && (i < nzz) && (j > 3) && (j < nxx-3) && (k > 3) && (k < nyy-3)) 
        {
            float dVx_dy = (75.0f*(Vx[i + j*nzz + (k-4)*nxx*nzz] - Vx[i + j*nzz + (k+3)*nxx*nzz]) +
                        1029.0f*(Vx[i + j*nzz + (k+2)*nxx*nzz] - Vx[i + j*nzz + (k-3)*nxx*nzz]) +
                        8575.0f*(Vx[i + j*nzz + (k-2)*nxx*nzz] - Vx[i + j*nzz + (k+1)*nxx*nzz]) +
                        128625.0f*(Vx[i + j*nzz + k*nxx*nzz]     - Vx[i + j*nzz + (k-1)*nxx*nzz])) / (107520.0f*dy);

            float dVy_dx = (75.0f*(Vy[i + (j-4)*nzz + k*nxx*nzz] - Vy[i + (j+3)*nzz + k*nxx*nzz]) +
                        1029.0f*(Vy[i + (j+2)*nzz + k*nxx*nzz] - Vy[i + (j-3)*nzz + k*nxx*nzz]) +
                        8575.0f*(Vy[i + (j-2)*nzz + k*nxx*nzz] - Vy[i + (j+1)*nzz + k*nxx*nzz]) +
                        128625.0f*(Vy[i + j*nzz + k*nxx*nzz]     - Vy[i + (j-1)*nzz + k*nxx*nzz])) / (107520.0f*dx);

            float C66xy = powf(0.25f*(1.0f/C66[i + (j+1)*nzz + (k+1)*nxx*nzz] + 1.0f/C66[i + (j+1)*nzz + k*nxx*nzz] + 
                                    1.0f/C66[i + j*nzz + (k+1)*nxx*nzz]     + 1.0f/C66[i + j*nzz + k*nxx*nzz]), -1.0f);

            Txy[index] += dt*C66xy*(dVx_dy + dVy_dx);
        }

        if((i > 3) && (i < nzz-3) && (j > 3) && (j < nxx-3) && (k >= 0) && (k < nyy)) 
        {
            float dVx_dz = (75.0f*(Vx[(i-4) + j*nzz + k*nxx*nzz] - Vx[(i+3) + j*nzz + k*nxx*nzz]) +
                        1029.0f*(Vx[(i+2) + j*nzz + k*nxx*nzz] - Vx[(i-3) + j*nzz + k*nxx*nzz]) +
                        8575.0f*(Vx[(i-2) + j*nzz + k*nxx*nzz] - Vx[(i+1) + j*nzz + k*nxx*nzz]) +
                        128625.0f*(Vx[i + j*nzz + k*nxx*nzz]     - Vx[(i-1) + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

            float dVz_dx = (75.0f*(Vz[i + (j-4)*nzz + k*nxx*nzz] - Vz[i + (j+3)*nzz + k*nxx*nzz]) +
                        1029.0f*(Vz[i + (j+2)*nzz + k*nxx*nzz] - Vz[i + (j-3)*nzz + k*nxx*nzz]) +
                        8575.0f*(Vz[i + (j-2)*nzz + k*nxx*nzz] - Vz[i + (j+1)*nzz + k*nxx*nzz]) +
                        128625.0f*(Vz[i + j*nzz + k*nxx*nzz]     - Vz[i + (j-1)*nzz + k*nxx*nzz])) / (107520.0f*dx);

            float C55xz = powf(0.25f*(1.0f/C55[(i+1) + (j+1)*nzz + k*nxx*nzz] + 1.0f/C55[i + (j+1)*nzz + k*nxx*nzz] + 
                                    1.0f/C55[(i+1) + j*nzz + k*nxx*nzz]     + 1.0f/C55[i + j*nzz + k*nxx*nzz]), -1.0f);

            Txz[index] += dt*C55xz*(dVx_dz + dVz_dx);
        }

        if((i > 3) && (i < nzz-3) && (j >= 0) && (j < nxx) && (k > 3) && (k < nyy-3)) 
        {
            float dVy_dz = (75.0f*(Vy[(i-4) + j*nzz + k*nxx*nzz] - Vy[(i+3) + j*nzz + k*nxx*nzz]) +
                        1029.0f*(Vy[(i+2) + j*nzz + k*nxx*nzz] - Vy[(i-3) + j*nzz + k*nxx*nzz]) +
                        8575.0f*(Vy[(i-2) + j*nzz + k*nxx*nzz] - Vy[(i+1) + j*nzz + k*nxx*nzz]) +
                        128625.0f*(Vy[i + j*nzz + k*nxx*nzz]     - Vy[(i-1) + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

            float dVz_dy = (75.0f*(Vz[i + j*nzz + (k-4)*nxx*nzz] - Vz[i + j*nzz + (k+3)*nxx*nzz]) +
                        1029.0f*(Vz[i + j*nzz + (k+2)*nxx*nzz] - Vz[i + j*nzz + (k-3)*nxx*nzz]) +
                        8575.0f*(Vz[i + j*nzz + (k-2)*nxx*nzz] - Vz[i + j*nzz + (k+1)*nxx*nzz]) +
                        128625.0f*(Vz[i + j*nzz + k*nxx*nzz]     - Vz[i + j*nzz + (k-1)*nxx*nzz])) / (107520.0f*dy);

            float C55yz = powf(0.25f*(1.0f/C55[(i+1) + j*nzz + (k+1)*nxx*nzz] + 1.0f/C55[i + j*nzz + (k+1)*nxx*nzz] + 
                                    1.0f/C55[(i+1) + j*nzz + k*nxx*nzz] +     1.0f/C55[i + j*nzz + k*nxx*nzz]), -1.0f);

            Tyz[index] += dt*C55yz*(dVy_dz + dVz_dy);
        }

        if ((i > 3) && (i < nzz-4) && (j > 3) && (j < nxx-4) && (k > 3) && (k < nyy-4))
        {
            P[index] = (Txx[index] + Tyy[index] + Tzz[index]) / 3.0f;
        }
    }
}

