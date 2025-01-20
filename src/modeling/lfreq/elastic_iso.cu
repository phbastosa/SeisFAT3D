# include "elastic_iso.cuh"

void Elastic_ISO::set_conditions()
{
    modeling_type = "elastic_iso";
    modeling_name = "Modeling type: Elastic isotropic solver";

    eikonal = new Eikonal_ISO();
    eikonal->parameters = parameters;
    eikonal->set_parameters();

    M = new float[volsize]();
    L = new float[volsize]();
    B = new float[volsize]();

    for (int index = 0; index < volsize; index++)
    {
        M[index] = Ro[index]*Vs[index]*Vs[index];
        L[index] = Ro[index]*Vp[index]*Vp[index] - 2.0f*M[index];
        B[index] = 1.0f / Ro[index];
    }
    
    cudaMalloc((void**)&(d_M), volsize*sizeof(float));
    cudaMalloc((void**)&(d_L), volsize*sizeof(float));
    cudaMalloc((void**)&(d_B), volsize*sizeof(float));

    cudaMemcpy(d_M, M, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_L, L, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, B, volsize*sizeof(float), cudaMemcpyHostToDevice);
}

void Elastic_ISO::propagation()
{
    for (int tId = 0; tId < nt + tlag; tId++)
    {
        compute_pressure<<<nBlocks, nThreads>>>(d_Vx, d_Vy, d_Vz, d_Txx, d_Tyy, d_Tzz, d_Txz, d_Tyz, d_Txy, d_P, d_M, d_L, d_T, wavelet, sIdx, sIdy, sIdz, tId, tlag, nt, dx, dy, dz, dt, nxx, nyy, nzz);
        cudaDeviceSynchronize();

        compute_velocity<<<nBlocks, nThreads>>>(d_Vx, d_Vy, d_Vz, d_Txx, d_Tyy, d_Tzz, d_Txz, d_Tyz, d_Txy, d_B, d_T, d1D, d2D, d3D, dx, dy, dz, dt, tId, tlag, nxx, nyy, nzz, nb);
        cudaDeviceSynchronize();

        compute_seismogram<<<sBlocks, nThreads>>>(d_P, rIdx, rIdy, rIdz, seismogram, geometry->spread[srcId], tId, tlag, nt, nxx, nzz);     
        cudaDeviceSynchronize();
    }
}

__global__ void compute_pressure(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * P, float * M, float * L, float * T, float * wavelet, int sIdx, int sIdy, int sIdz, int tId, int tlag, int nt, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz)
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

            Txx[index] += dt*((L[index] + 2*M[index])*dVx_dx + L[index]*(dVy_dy + dVz_dz));
            Tyy[index] += dt*((L[index] + 2*M[index])*dVy_dy + L[index]*(dVx_dx + dVz_dz));
            Tzz[index] += dt*((L[index] + 2*M[index])*dVz_dz + L[index]*(dVx_dx + dVy_dy));                    
        }

        if((i >= 3) && (i < nzz-4) && (j > 3) && (j < nxx-3) && (k > 3) && (k < nyy-3)) 
        {
            float dVx_dy = (75.0f*(Vx[i + j*nzz + (k-4)*nxx*nzz] - Vx[i + j*nzz + (k+3)*nxx*nzz]) +
                          1029.0f*(Vx[i + j*nzz + (k+2)*nxx*nzz] - Vx[i + j*nzz + (k-3)*nxx*nzz]) +
                          8575.0f*(Vx[i + j*nzz + (k-2)*nxx*nzz] - Vx[i + j*nzz + (k+1)*nxx*nzz]) +
                        128625.0f*(Vx[i + j*nzz + k*nxx*nzz]     - Vx[i + j*nzz + (k-1)*nxx*nzz])) / (107520.0f*dy);

            float dVy_dx = (75.0f*(Vy[i + (j-4)*nzz + k*nxx*nzz] - Vy[i + (j+3)*nzz + k*nxx*nzz]) +
                          1029.0f*(Vy[i + (j+2)*nzz + k*nxx*nzz] - Vy[i + (j-3)*nzz + k*nxx*nzz]) +
                          8575.0f*(Vy[i + (j-2)*nzz + k*nxx*nzz] - Vy[i + (j+1)*nzz + k*nxx*nzz]) +
                        128625.0f*(Vy[i + j*nzz + k*nxx*nzz]     - Vy[i + (j-1)*nzz + k*nxx*nzz])) / (107520.0f*dx);

            float Mxy = powf(0.25f*(1.0f/M[i + (j+1)*nzz + (k+1)*nxx*nzz] + 1.0f/M[i + (j+1)*nzz + k*nxx*nzz] + 
                                    1.0f/M[i + j*nzz + (k+1)*nxx*nzz]     + 1.0f/M[i + j*nzz + k*nxx*nzz]), -1.0f);

            Txy[index] += dt*Mxy*(dVx_dy + dVy_dx);
        }

        if((i > 3) && (i < nzz-3) && (j > 3) && (j < nxx-3) && (k >= 3) && (k < nyy-4)) 
        {
            float dVx_dz = (75.0f*(Vx[(i-4) + j*nzz + k*nxx*nzz] - Vx[(i+3) + j*nzz + k*nxx*nzz]) +
                          1029.0f*(Vx[(i+2) + j*nzz + k*nxx*nzz] - Vx[(i-3) + j*nzz + k*nxx*nzz]) +
                          8575.0f*(Vx[(i-2) + j*nzz + k*nxx*nzz] - Vx[(i+1) + j*nzz + k*nxx*nzz]) +
                        128625.0f*(Vx[i + j*nzz + k*nxx*nzz]     - Vx[(i-1) + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

            float dVz_dx = (75.0f*(Vz[i + (j-4)*nzz + k*nxx*nzz] - Vz[i + (j+3)*nzz + k*nxx*nzz]) +
                          1029.0f*(Vz[i + (j+2)*nzz + k*nxx*nzz] - Vz[i + (j-3)*nzz + k*nxx*nzz]) +
                          8575.0f*(Vz[i + (j-2)*nzz + k*nxx*nzz] - Vz[i + (j+1)*nzz + k*nxx*nzz]) +
                        128625.0f*(Vz[i + j*nzz + k*nxx*nzz]     - Vz[i + (j-1)*nzz + k*nxx*nzz])) / (107520.0f*dx);

            float Mxz = powf(0.25f*(1.0f/M[(i+1) + (j+1)*nzz + k*nxx*nzz] + 1.0f/M[i + (j+1)*nzz + k*nxx*nzz] + 
                                    1.0f/M[(i+1) + j*nzz + k*nxx*nzz]     + 1.0f/M[i + j*nzz + k*nxx*nzz]), -1.0f);

            Txz[index] += dt*Mxz*(dVx_dz + dVz_dx);
        }

        if((i > 3) && (i < nzz-3) && (j >= 3) && (j < nxx-4) && (k > 3) && (k < nyy-3)) 
        {
            float dVy_dz = (75.0f*(Vy[(i-4) + j*nzz + k*nxx*nzz] - Vy[(i+3) + j*nzz + k*nxx*nzz]) +
                          1029.0f*(Vy[(i+2) + j*nzz + k*nxx*nzz] - Vy[(i-3) + j*nzz + k*nxx*nzz]) +
                          8575.0f*(Vy[(i-2) + j*nzz + k*nxx*nzz] - Vy[(i+1) + j*nzz + k*nxx*nzz]) +
                        128625.0f*(Vy[i + j*nzz + k*nxx*nzz]     - Vy[(i-1) + j*nzz + k*nxx*nzz])) / (107520.0f*dz);

            float dVz_dy = (75.0f*(Vz[i + j*nzz + (k-4)*nxx*nzz] - Vz[i + j*nzz + (k+3)*nxx*nzz]) +
                          1029.0f*(Vz[i + j*nzz + (k+2)*nxx*nzz] - Vz[i + j*nzz + (k-3)*nxx*nzz]) +
                          8575.0f*(Vz[i + j*nzz + (k-2)*nxx*nzz] - Vz[i + j*nzz + (k+1)*nxx*nzz]) +
                        128625.0f*(Vz[i + j*nzz + k*nxx*nzz]     - Vz[i + j*nzz + (k-1)*nxx*nzz])) / (107520.0f*dy);

            float Myz = powf(0.25f*(1.0f/M[(i+1) + j*nzz + (k+1)*nxx*nzz] + 1.0f/M[i + j*nzz + (k+1)*nxx*nzz] + 
                                    1.0f/M[(i+1) + j*nzz + k*nxx*nzz] +     1.0f/M[i + j*nzz + k*nxx*nzz]), -1.0f);

            Tyz[index] += dt*Myz*(dVy_dz + dVz_dy);
        }

        if ((i > 3) && (i < nzz-4) && (j > 3) && (j < nxx-4) && (k > 3) && (k < nyy-4))
        {
            P[index] = (Txx[index] + Tyy[index] + Tzz[index]) / 3.0f;
        }
    }
}

