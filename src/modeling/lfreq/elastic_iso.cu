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
        compute_velocity_ssg<<<nBlocks, nThreads>>>(d_Vx, d_Vy, d_Vz, d_Txx, d_Tyy, d_Tzz, d_Txz, d_Tyz, d_Txy, d_B, d_T, d1D, d2D, d3D, wavelet, dx, dy, dz, dt, tId, tlag, sIdx, sIdy, sIdz, nxx, nyy, nzz, nb, nt);

        compute_pressure_ssg<<<nBlocks, nThreads>>>(d_Vx, d_Vy, d_Vz, d_Txx, d_Tyy, d_Tzz, d_Txz, d_Tyz, d_Txy, d_P, d_M, d_L, d_T, tId, tlag, dx, dy, dz, dt, nxx, nyy, nzz);

        compute_seismogram<<<sBlocks, nThreads>>>(d_P, rIdx, rIdy, rIdz, seismogram, geometry->spread[srcId], tId, tlag, nt, nxx, nzz);     
    }
}

__global__ void compute_velocity_ssg(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * B, float * T,  
                                     float * damp1D, float * damp2D, float * damp3D, float * wavelet, float dx, float dy, float dz, float dt, int tId, int tlag, int sIdx, 
                                     int sIdy, int sIdz, int nxx, int nyy, int nzz, int nb, int nt)
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

    if ((T[index] < (float)(tId + tlag)*dt) && (index < nxx*nyy*nzz))
    {
        if((i >= 3) && (i < nzz-4) && (j > 3) && (j < nxx-3) && (k >= 3) && (k < nyy-4)) 
        {
            float dTxx_dx = (FDM1*(Txx[i + (j-4)*nzz + k*nxx*nzz] - Txx[i + (j+3)*nzz + k*nxx*nzz]) +
                             FDM2*(Txx[i + (j+2)*nzz + k*nxx*nzz] - Txx[i + (j-3)*nzz + k*nxx*nzz]) +
                             FDM3*(Txx[i + (j-2)*nzz + k*nxx*nzz] - Txx[i + (j+1)*nzz + k*nxx*nzz]) +
                             FDM4*(Txx[i + j*nzz + k*nxx*nzz]     - Txx[i + (j-1)*nzz + k*nxx*nzz])) / dx;

            float dTxy_dy = (FDM1*(Txy[i + j*nzz + (k-3)*nxx*nzz] - Txy[i + j*nzz + (k+4)*nxx*nzz]) +
                             FDM2*(Txy[i + j*nzz + (k+3)*nxx*nzz] - Txy[i + j*nzz + (k-2)*nxx*nzz]) +
                             FDM3*(Txy[i + j*nzz + (k-1)*nxx*nzz] - Txy[i + j*nzz + (k+2)*nxx*nzz]) +
                             FDM4*(Txy[i + j*nzz + (k+1)*nxx*nzz] - Txy[i + j*nzz + k*nxx*nzz])) / dy;

            float dTxz_dz = (FDM1*(Txz[(i-3) + j*nzz + k*nxx*nzz] - Txz[(i+4) + j*nzz + k*nxx*nzz]) +
                             FDM2*(Txz[(i+3) + j*nzz + k*nxx*nzz] - Txz[(i-2) + j*nzz + k*nxx*nzz]) +
                             FDM3*(Txz[(i-1) + j*nzz + k*nxx*nzz] - Txz[(i+2) + j*nzz + k*nxx*nzz]) +
                             FDM4*(Txz[(i+1) + j*nzz + k*nxx*nzz] - Txz[i + j*nzz + k*nxx*nzz])) / dz;

            float Bx = 0.5f*(B[i + (j+1)*nzz + k*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);

            Vx[index] += dt*Bx*(dTxx_dx + dTxy_dy + dTxz_dz); 
        }

        if((i >= 3) && (i < nzz-3) && (j >= 3) && (j < nxx-4) && (k > 3) && (k < nyy-3)) 
        {
            float dTxy_dx = (FDM1*(Txy[i + (j-3)*nzz + k*nxx*nzz] - Txy[i + (j+4)*nzz + k*nxx*nzz]) +
                             FDM2*(Txy[i + (j+3)*nzz + k*nxx*nzz] - Txy[i + (j-2)*nzz + k*nxx*nzz]) +
                             FDM3*(Txy[i + (j-1)*nzz + k*nxx*nzz] - Txy[i + (j+2)*nzz + k*nxx*nzz]) +
                             FDM4*(Txy[i + (j+1)*nzz + k*nxx*nzz] - Txy[i + j*nzz + k*nxx*nzz])) / dx;

            float dTyy_dy = (FDM1*(Tyy[i + j*nzz + (k-4)*nxx*nzz] - Tyy[i + j*nzz + (k+3)*nxx*nzz]) +
                             FDM2*(Tyy[i + j*nzz + (k+2)*nxx*nzz] - Tyy[i + j*nzz + (k-3)*nxx*nzz]) +
                             FDM3*(Tyy[i + j*nzz + (k-2)*nxx*nzz] - Tyy[i + j*nzz + (k+1)*nxx*nzz]) +
                             FDM4*(Tyy[i + j*nzz + k*nxx*nzz]     - Tyy[i + j*nzz + (k-1)*nxx*nzz])) / dy;

            float dTyz_dz = (FDM1*(Tyz[(i-3) + j*nzz + k*nxx*nzz] - Tyz[(i+4) + j*nzz + k*nxx*nzz]) +
                             FDM2*(Tyz[(i+3) + j*nzz + k*nxx*nzz] - Tyz[(i-2) + j*nzz + k*nxx*nzz]) +
                             FDM3*(Tyz[(i-1) + j*nzz + k*nxx*nzz] - Tyz[(i+2) + j*nzz + k*nxx*nzz]) +
                             FDM4*(Tyz[(i+1) + j*nzz + k*nxx*nzz] - Tyz[i + j*nzz + k*nxx*nzz])) / dz;

            float By = 0.5f*(B[i + j*nzz + (k+1)*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);

            Vy[index] += dt*By*(dTxy_dx + dTyy_dy + dTyz_dz); 
        }    

        if((i > 3) && (i < nzz-3) && (j >= 3) && (j < nxx-4) && (k >= 3) && (k < nyy-4)) 
        {
            float dTxz_dx = (FDM1*(Txz[i + (j-3)*nzz + k*nxx*nzz] - Txz[i + (j+4)*nzz + k*nxx*nzz]) +
                             FDM2*(Txz[i + (j+3)*nzz + k*nxx*nzz] - Txz[i + (j-2)*nzz + k*nxx*nzz]) +
                             FDM3*(Txz[i + (j-1)*nzz + k*nxx*nzz] - Txz[i + (j+2)*nzz + k*nxx*nzz]) +
                             FDM4*(Txz[i + (j+1)*nzz + k*nxx*nzz] - Txz[i + j*nzz + k*nxx*nzz])) / dx;

            float dTyz_dy = (FDM1*(Tyz[i + j*nzz + (k-3)*nxx*nzz] - Tyz[i + j*nzz + (k+4)*nxx*nzz]) +
                             FDM2*(Tyz[i + j*nzz + (k+3)*nxx*nzz] - Tyz[i + j*nzz + (k-2)*nxx*nzz]) +
                             FDM3*(Tyz[i + j*nzz + (k-1)*nxx*nzz] - Tyz[i + j*nzz + (k+2)*nxx*nzz]) +
                             FDM4*(Tyz[i + j*nzz + (k+1)*nxx*nzz] - Tyz[i + j*nzz + k*nxx*nzz])) / dy;

            float dTzz_dz = (FDM1*(Tzz[(i-4) + j*nzz + k*nxx*nzz] - Tzz[(i+3) + j*nzz + k*nxx*nzz]) +
                             FDM2*(Tzz[(i+2) + j*nzz + k*nxx*nzz] - Tzz[(i-3) + j*nzz + k*nxx*nzz]) +
                             FDM3*(Tzz[(i-2) + j*nzz + k*nxx*nzz] - Tzz[(i+1) + j*nzz + k*nxx*nzz]) +
                             FDM4*(Tzz[i + j*nzz + k*nxx*nzz]     - Tzz[(i-1) + j*nzz + k*nxx*nzz])) / dz;

            float Bz = 0.5f*(B[(i+1) + j*nzz + k*nxx*nzz] + B[i + j*nzz + k*nxx*nzz]);

            Vz[index] += dt*Bz*(dTxz_dx + dTyz_dy + dTzz_dz); 
        }

    	float damper = get_boundary_damper(damp1D, damp2D, damp3D, i, j, k, nxx, nyy, nzz, nb);

        Vx[index] *= damper;
        Vy[index] *= damper;
        Vz[index] *= damper;

        Txx[index] *= damper;
        Tyy[index] *= damper;
        Tzz[index] *= damper;
        Txz[index] *= damper;
        Tyz[index] *= damper;
        Txy[index] *= damper;
    }
}

__global__ void compute_pressure_ssg(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * P, float * M, 
                                     float * L, float * T, int tId, int tlag, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         
    int j = (int) (index - k*nxx*nzz) / nzz;   
    int i = (int) (index - j*nzz - k*nxx*nzz); 

    if ((T[index] < (float)(tId + tlag)*dt) && (index < nxx*nyy*nzz))
    {
        if((i >= 3) && (i < nzz-4) && (j >= 3) && (j < nxx-4) && (k >= 3) && (k < nyy-4)) 
        {    
            float dVx_dx = (FDM1*(Vx[i + (j-3)*nzz + k*nxx*nzz] - Vx[i + (j+4)*nzz + k*nxx*nzz]) +
                            FDM2*(Vx[i + (j+3)*nzz + k*nxx*nzz] - Vx[i + (j-2)*nzz + k*nxx*nzz]) +
                            FDM3*(Vx[i + (j-1)*nzz + k*nxx*nzz] - Vx[i + (j+2)*nzz + k*nxx*nzz]) +
                            FDM4*(Vx[i + (j+1)*nzz + k*nxx*nzz] - Vx[i + j*nzz + k*nxx*nzz])) / dx;

            float dVy_dy = (FDM1*(Vy[i + j*nzz + (k-3)*nxx*nzz] - Vy[i + j*nzz + (k+4)*nxx*nzz]) +
                            FDM2*(Vy[i + j*nzz + (k+3)*nxx*nzz] - Vy[i + j*nzz + (k-2)*nxx*nzz]) +
                            FDM3*(Vy[i + j*nzz + (k-1)*nxx*nzz] - Vy[i + j*nzz + (k+2)*nxx*nzz]) +
                            FDM4*(Vy[i + j*nzz + (k+1)*nxx*nzz] - Vy[i + j*nzz + k*nxx*nzz])) / dy;

            float dVz_dz = (FDM1*(Vz[(i-3) + j*nzz + k*nxx*nzz] - Vz[(i+4) + j*nzz + k*nxx*nzz]) +
                            FDM2*(Vz[(i+3) + j*nzz + k*nxx*nzz] - Vz[(i-2) + j*nzz + k*nxx*nzz]) +
                            FDM3*(Vz[(i-1) + j*nzz + k*nxx*nzz] - Vz[(i+2) + j*nzz + k*nxx*nzz]) +
                            FDM4*(Vz[(i+1) + j*nzz + k*nxx*nzz] - Vz[i + j*nzz + k*nxx*nzz])) / dz;

            Txx[index] += dt*((L[index] + 2*M[index])*dVx_dx + L[index]*(dVy_dy + dVz_dz));
            Tyy[index] += dt*((L[index] + 2*M[index])*dVy_dy + L[index]*(dVx_dx + dVz_dz));
            Tzz[index] += dt*((L[index] + 2*M[index])*dVz_dz + L[index]*(dVx_dx + dVy_dy));                    
        }

        if((i >= 3) && (i < nzz-4) && (j > 3) && (j < nxx-3) && (k > 3) && (k < nyy-3)) 
        {
            float dVx_dy = (FDM1*(Vx[i + j*nzz + (k-4)*nxx*nzz] - Vx[i + j*nzz + (k+3)*nxx*nzz]) +
                            FDM2*(Vx[i + j*nzz + (k+2)*nxx*nzz] - Vx[i + j*nzz + (k-3)*nxx*nzz]) +
                            FDM3*(Vx[i + j*nzz + (k-2)*nxx*nzz] - Vx[i + j*nzz + (k+1)*nxx*nzz]) +
                            FDM4*(Vx[i + j*nzz + k*nxx*nzz]     - Vx[i + j*nzz + (k-1)*nxx*nzz])) / dy;

            float dVy_dx = (FDM1*(Vy[i + (j-4)*nzz + k*nxx*nzz] - Vy[i + (j+3)*nzz + k*nxx*nzz]) +
                            FDM2*(Vy[i + (j+2)*nzz + k*nxx*nzz] - Vy[i + (j-3)*nzz + k*nxx*nzz]) +
                            FDM3*(Vy[i + (j-2)*nzz + k*nxx*nzz] - Vy[i + (j+1)*nzz + k*nxx*nzz]) +
                            FDM4*(Vy[i + j*nzz + k*nxx*nzz]     - Vy[i + (j-1)*nzz + k*nxx*nzz])) / dx;

            float Mxy = powf(0.25f*(1.0f/M[i + (j+1)*nzz + (k+1)*nxx*nzz] + 1.0f/M[i + (j+1)*nzz + k*nxx*nzz] + 
                                    1.0f/M[i + j*nzz + (k+1)*nxx*nzz]     + 1.0f/M[i + j*nzz + k*nxx*nzz]), -1.0f);

            Txy[index] += dt*Mxy*(dVx_dy + dVy_dx);
        }

        if((i > 3) && (i < nzz-3) && (j > 3) && (j < nxx-3) && (k >= 3) && (k < nyy-4)) 
        {
            float dVx_dz = (FDM1*(Vx[(i-4) + j*nzz + k*nxx*nzz] - Vx[(i+3) + j*nzz + k*nxx*nzz]) +
                            FDM2*(Vx[(i+2) + j*nzz + k*nxx*nzz] - Vx[(i-3) + j*nzz + k*nxx*nzz]) +
                            FDM3*(Vx[(i-2) + j*nzz + k*nxx*nzz] - Vx[(i+1) + j*nzz + k*nxx*nzz]) +
                            FDM4*(Vx[i + j*nzz + k*nxx*nzz]     - Vx[(i-1) + j*nzz + k*nxx*nzz])) / dz;

            float dVz_dx = (FDM1*(Vz[i + (j-4)*nzz + k*nxx*nzz] - Vz[i + (j+3)*nzz + k*nxx*nzz]) +
                            FDM2*(Vz[i + (j+2)*nzz + k*nxx*nzz] - Vz[i + (j-3)*nzz + k*nxx*nzz]) +
                            FDM3*(Vz[i + (j-2)*nzz + k*nxx*nzz] - Vz[i + (j+1)*nzz + k*nxx*nzz]) +
                            FDM4*(Vz[i + j*nzz + k*nxx*nzz]     - Vz[i + (j-1)*nzz + k*nxx*nzz])) / dx;

            float Mxz = powf(0.25f*(1.0f/M[(i+1) + (j+1)*nzz + k*nxx*nzz] + 1.0f/M[i + (j+1)*nzz + k*nxx*nzz] + 
                                    1.0f/M[(i+1) + j*nzz + k*nxx*nzz]     + 1.0f/M[i + j*nzz + k*nxx*nzz]), -1.0f);

            Txz[index] += dt*Mxz*(dVx_dz + dVz_dx);
        }

        if((i > 3) && (i < nzz-3) && (j >= 3) && (j < nxx-4) && (k > 3) && (k < nyy-3)) 
        {
            float dVy_dz = (FDM1*(Vy[(i-4) + j*nzz + k*nxx*nzz] - Vy[(i+3) + j*nzz + k*nxx*nzz]) +
                            FDM2*(Vy[(i+2) + j*nzz + k*nxx*nzz] - Vy[(i-3) + j*nzz + k*nxx*nzz]) +
                            FDM3*(Vy[(i-2) + j*nzz + k*nxx*nzz] - Vy[(i+1) + j*nzz + k*nxx*nzz]) +
                            FDM4*(Vy[i + j*nzz + k*nxx*nzz]     - Vy[(i-1) + j*nzz + k*nxx*nzz])) / dz;

            float dVz_dy = (FDM1*(Vz[i + j*nzz + (k-4)*nxx*nzz] - Vz[i + j*nzz + (k+3)*nxx*nzz]) +
                            FDM2*(Vz[i + j*nzz + (k+2)*nxx*nzz] - Vz[i + j*nzz + (k-3)*nxx*nzz]) +
                            FDM3*(Vz[i + j*nzz + (k-2)*nxx*nzz] - Vz[i + j*nzz + (k+1)*nxx*nzz]) +
                            FDM4*(Vz[i + j*nzz + k*nxx*nzz]     - Vz[i + j*nzz + (k-1)*nxx*nzz])) / dy;

            float Myz = powf(0.25f*(1.0f/M[(i+1) + j*nzz + (k+1)*nxx*nzz] + 1.0f/M[i + j*nzz + (k+1)*nxx*nzz] + 
                                    1.0f/M[(i+1) + j*nzz + k*nxx*nzz] +     1.0f/M[i + j*nzz + k*nxx*nzz]), -1.0f);

            Tyz[index] += dt*Myz*(dVy_dz + dVz_dy);
        }

        if ((i > 3) && (i < nzz-4) && (j > 3) && (j < nxx-4) && (k > 3) && (k < nyy-4))
            P[index] = (Txx[index] + Tyy[index] + Tzz[index]) / 3.0f;
    }
}

