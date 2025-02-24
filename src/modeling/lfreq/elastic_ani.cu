# include "elastic_ani.cuh"

void Elastic_ANI::set_conditions()
{
    modeling_type = "elastic_ani";
    modeling_name = "Modeling type: Elastic anisotropic solver";

    eikonal = new Eikonal_ANI();
    eikonal->parameters = parameters;
    eikonal->set_parameters();

    B = new float[volsize]();

    # pragma omp parallel for
    for (int index = 0; index < volsize; index++)
        B[index] = 1.0f / Ro[index];

    cudaMalloc((void**)&(d_B), volsize*sizeof(float));
    cudaMemcpy(d_B, B, volsize*sizeof(float), cudaMemcpyHostToDevice);

    float * Cij = new float[nPoints]();

    std::string Cijkl_folder = catch_parameter("Cijkl_folder", parameters);

    float * C11 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C11.bin", Cij, nPoints);
    expand_boundary(Cij, C11);
    cudaMalloc((void**)&(d_C11), volsize*sizeof(float));
    cudaMemcpy(d_C11, C11, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C11;

    float * C12 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C12.bin", Cij, nPoints);
    expand_boundary(Cij, C12);
    cudaMalloc((void**)&(d_C12), volsize*sizeof(float));
    cudaMemcpy(d_C12, C12, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C12;

    float * C13 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C13.bin", Cij, nPoints);
    expand_boundary(Cij, C13);
    cudaMalloc((void**)&(d_C13), volsize*sizeof(float));
    cudaMemcpy(d_C13, C13, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C13;

    float * C14 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C14.bin", Cij, nPoints);
    expand_boundary(Cij, C14);
    cudaMalloc((void**)&(d_C14), volsize*sizeof(float));
    cudaMemcpy(d_C14, C14, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C14;

    float * C15 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C15.bin", Cij, nPoints);
    expand_boundary(Cij, C15);
    cudaMalloc((void**)&(d_C15), volsize*sizeof(float));
    cudaMemcpy(d_C15, C15, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C15;

    float * C16 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C16.bin", Cij, nPoints);
    expand_boundary(Cij, C16);
    cudaMalloc((void**)&(d_C16), volsize*sizeof(float));
    cudaMemcpy(d_C16, C16, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C16;

    float * C22 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C22.bin", Cij, nPoints);
    expand_boundary(Cij, C22);
    cudaMalloc((void**)&(d_C22), volsize*sizeof(float));
    cudaMemcpy(d_C22, C22, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C22;

    float * C23 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C23.bin", Cij, nPoints);
    expand_boundary(Cij, C23);
    cudaMalloc((void**)&(d_C23), volsize*sizeof(float));
    cudaMemcpy(d_C23, C23, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C23;
    
    float * C24 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C24.bin", Cij, nPoints);
    expand_boundary(Cij, C24);
    cudaMalloc((void**)&(d_C24), volsize*sizeof(float));
    cudaMemcpy(d_C24, C24, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C24;

    float * C25 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C25.bin", Cij, nPoints);
    expand_boundary(Cij, C25);
    cudaMalloc((void**)&(d_C25), volsize*sizeof(float));
    cudaMemcpy(d_C25, C25, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C25;

    float * C26 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C26.bin", Cij, nPoints);
    expand_boundary(Cij, C26);
    cudaMalloc((void**)&(d_C26), volsize*sizeof(float));
    cudaMemcpy(d_C26, C26, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C26;

    float * C33 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C33.bin", Cij, nPoints);
    expand_boundary(Cij, C33);
    cudaMalloc((void**)&(d_C33), volsize*sizeof(float));
    cudaMemcpy(d_C33, C33, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C33;
    
    float * C34 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C34.bin", Cij, nPoints);
    expand_boundary(Cij, C34);
    cudaMalloc((void**)&(d_C34), volsize*sizeof(float));
    cudaMemcpy(d_C34, C34, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C34;

    float * C35 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C35.bin", Cij, nPoints);
    expand_boundary(Cij, C35);
    cudaMalloc((void**)&(d_C35), volsize*sizeof(float));
    cudaMemcpy(d_C35, C35, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C35;

    float * C36 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C36.bin", Cij, nPoints);
    expand_boundary(Cij, C36);
    cudaMalloc((void**)&(d_C36), volsize*sizeof(float));
    cudaMemcpy(d_C36, C36, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C36;

    float * C44 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C44.bin", Cij, nPoints);
    expand_boundary(Cij, C44);
    cudaMalloc((void**)&(d_C44), volsize*sizeof(float));
    cudaMemcpy(d_C44, C44, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C44;

    float * C45 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C45.bin", Cij, nPoints);
    expand_boundary(Cij, C45);
    cudaMalloc((void**)&(d_C45), volsize*sizeof(float));
    cudaMemcpy(d_C45, C45, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C45;

    float * C46 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C46.bin", Cij, nPoints);
    expand_boundary(Cij, C46);
    cudaMalloc((void**)&(d_C46), volsize*sizeof(float));
    cudaMemcpy(d_C46, C46, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C46;

    float * C55 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C55.bin", Cij, nPoints);
    expand_boundary(Cij, C55);
    cudaMalloc((void**)&(d_C55), volsize*sizeof(float));
    cudaMemcpy(d_C55, C55, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C55;

    float * C56 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C56.bin", Cij, nPoints);
    expand_boundary(Cij, C56);
    cudaMalloc((void**)&(d_C56), volsize*sizeof(float));
    cudaMemcpy(d_C56, C56, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C56;

    float * C66 = new float[volsize]();
    import_binary_float(Cijkl_folder + "C66.bin", Cij, nPoints);
    expand_boundary(Cij, C66);
    cudaMalloc((void**)&(d_C66), volsize*sizeof(float));
    cudaMemcpy(d_C66, C66, volsize*sizeof(float), cudaMemcpyHostToDevice);
    delete[] C66;
    
    delete[] Cij;
}

void Elastic_ANI::propagation()
{
    for (int tId = 0; tId < nt + tlag; tId++)
    {
        compute_pressure<<<nBlocks, nThreads>>>(d_Vx, d_Vy, d_Vz, d_Txx, d_Tyy, d_Tzz, d_Txz, d_Tyz, d_Txy, d_P, d_T, 
                                                d_C11, d_C12, d_C13, d_C14, d_C15, d_C16, d_C22, d_C23, d_C24, d_C25, d_C26, 
                                                d_C33, d_C34, d_C35, d_C36, d_C44, d_C45, d_C46, d_C55, d_C56, d_C66, wavelet, 
                                                sIdx, sIdy, sIdz, tId, tlag, nt, dx, dy, dz, dt, nxx, nyy, nzz);
        cudaDeviceSynchronize();

        compute_velocity<<<nBlocks, nThreads>>>(d_Vx, d_Vy, d_Vz, d_Txx, d_Tyy, d_Tzz, d_Txz, d_Tyz, d_Txy, d_B, d_T, d1D, d2D, d3D, dx, dy, dz, dt, tId, tlag, nxx, nyy, nzz, nb);
        cudaDeviceSynchronize();

        compute_seismogram<<<sBlocks, nThreads>>>(d_P, rIdx, rIdy, rIdz, seismogram, geometry->spread[srcId], tId, tlag, nt, nxx, nzz);     
        cudaDeviceSynchronize();
    }
}

__global__ void compute_pressure(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * P, float * T, 
                                 float * C11, float * C12, float * C13, float * C14, float * C15, float * C16, float * C22, float * C23, float * C24, float * C25, float * C26,  
                                 float * C33, float * C34, float * C35, float * C36, float * C44, float * C45, float * C46, float * C55, float * C56, float * C66, float * wavelet, 
                                 int sIdx, int sIdy, int sIdz, int tId, int tlag, int nt, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz)
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
    // if ((index < nxx*nyy*nzz))
    {
        if((i >= 3) && (i < nzz-4) && (j >= 3) && (j < nxx-4) && (k >= 3) && (k < nyy-4)) 
        {            

            float dVx_dx = (FDM1*(Vx[i + (j-3)*nzz + k*nxx*nzz] - Vx[i + (j+4)*nzz + k*nxx*nzz]) +
                            FDM2*(Vx[i + (j+3)*nzz + k*nxx*nzz] - Vx[i + (j-2)*nzz + k*nxx*nzz]) +
                            FDM3*(Vx[i + (j-1)*nzz + k*nxx*nzz] - Vx[i + (j+2)*nzz + k*nxx*nzz]) +
                            FDM4*(Vx[i + (j+1)*nzz + k*nxx*nzz] - Vx[i + j*nzz + k*nxx*nzz])) / dx;

            float dVy_dx = (FDM1*(Vy[i + (j-3)*nzz + k*nxx*nzz] - Vy[i + (j+4)*nzz + k*nxx*nzz]) +
                            FDM2*(Vy[i + (j+3)*nzz + k*nxx*nzz] - Vy[i + (j-2)*nzz + k*nxx*nzz]) +
                            FDM3*(Vy[i + (j-1)*nzz + k*nxx*nzz] - Vy[i + (j+2)*nzz + k*nxx*nzz]) +
                            FDM4*(Vy[i + (j+1)*nzz + k*nxx*nzz] - Vy[i + j*nzz + k*nxx*nzz])) / dx;

            float dVz_dx = (FDM1*(Vz[i + (j-3)*nzz + k*nxx*nzz] - Vz[i + (j+4)*nzz + k*nxx*nzz]) +
                            FDM2*(Vz[i + (j+3)*nzz + k*nxx*nzz] - Vz[i + (j-2)*nzz + k*nxx*nzz]) +
                            FDM3*(Vz[i + (j-1)*nzz + k*nxx*nzz] - Vz[i + (j+2)*nzz + k*nxx*nzz]) +
                            FDM4*(Vz[i + (j+1)*nzz + k*nxx*nzz] - Vz[i + j*nzz + k*nxx*nzz])) / dx;

            float dVx_dy = (FDM1*(Vx[i + j*nzz + (k-3)*nxx*nzz] - Vx[i + j*nzz + (k+4)*nxx*nzz]) +
                            FDM2*(Vx[i + j*nzz + (k+3)*nxx*nzz] - Vx[i + j*nzz + (k-2)*nxx*nzz]) +
                            FDM3*(Vx[i + j*nzz + (k-1)*nxx*nzz] - Vx[i + j*nzz + (k+2)*nxx*nzz]) +
                            FDM4*(Vx[i + j*nzz + (k+1)*nxx*nzz] - Vx[i + j*nzz + k*nxx*nzz])) / dy;

            float dVy_dy = (FDM1*(Vy[i + j*nzz + (k-3)*nxx*nzz] - Vy[i + j*nzz + (k+4)*nxx*nzz]) +
                            FDM2*(Vy[i + j*nzz + (k+3)*nxx*nzz] - Vy[i + j*nzz + (k-2)*nxx*nzz]) +
                            FDM3*(Vy[i + j*nzz + (k-1)*nxx*nzz] - Vy[i + j*nzz + (k+2)*nxx*nzz]) +
                            FDM4*(Vy[i + j*nzz + (k+1)*nxx*nzz] - Vy[i + j*nzz + k*nxx*nzz])) / dy;

            float dVz_dy = (FDM1*(Vz[i + j*nzz + (k-3)*nxx*nzz] - Vz[i + j*nzz + (k+4)*nxx*nzz]) +
                            FDM2*(Vz[i + j*nzz + (k+3)*nxx*nzz] - Vz[i + j*nzz + (k-2)*nxx*nzz]) +
                            FDM3*(Vz[i + j*nzz + (k-1)*nxx*nzz] - Vz[i + j*nzz + (k+2)*nxx*nzz]) +
                            FDM4*(Vz[i + j*nzz + (k+1)*nxx*nzz] - Vz[i + j*nzz + k*nxx*nzz])) / dy;

            float dVx_dz = (FDM1*(Vx[(i-3) + j*nzz + k*nxx*nzz] - Vx[(i+4) + j*nzz + k*nxx*nzz]) +
                            FDM2*(Vx[(i+3) + j*nzz + k*nxx*nzz] - Vx[(i-2) + j*nzz + k*nxx*nzz]) +
                            FDM3*(Vx[(i-1) + j*nzz + k*nxx*nzz] - Vx[(i+2) + j*nzz + k*nxx*nzz]) +
                            FDM4*(Vx[(i+1) + j*nzz + k*nxx*nzz] - Vx[i + j*nzz + k*nxx*nzz])) / dz;

            float dVy_dz = (FDM1*(Vy[(i-3) + j*nzz + k*nxx*nzz] - Vy[(i+4) + j*nzz + k*nxx*nzz]) +
                            FDM2*(Vy[(i+3) + j*nzz + k*nxx*nzz] - Vy[(i-2) + j*nzz + k*nxx*nzz]) +
                            FDM3*(Vy[(i-1) + j*nzz + k*nxx*nzz] - Vy[(i+2) + j*nzz + k*nxx*nzz]) +
                            FDM4*(Vy[(i+1) + j*nzz + k*nxx*nzz] - Vy[i + j*nzz + k*nxx*nzz])) / dz;

            float dVz_dz = (FDM1*(Vz[(i-3) + j*nzz + k*nxx*nzz] - Vz[(i+4) + j*nzz + k*nxx*nzz]) +
                            FDM2*(Vz[(i+3) + j*nzz + k*nxx*nzz] - Vz[(i-2) + j*nzz + k*nxx*nzz]) +
                            FDM3*(Vz[(i-1) + j*nzz + k*nxx*nzz] - Vz[(i+2) + j*nzz + k*nxx*nzz]) +
                            FDM4*(Vz[(i+1) + j*nzz + k*nxx*nzz] - Vz[i + j*nzz + k*nxx*nzz])) / dz;

            Txx[index] += dt*(C11[index]*dVx_dx + C16[index]*dVx_dy + C15[index]*dVx_dz +
                              C16[index]*dVy_dx + C12[index]*dVy_dy + C14[index]*dVy_dz +
                              C15[index]*dVz_dx + C14[index]*dVz_dy + C13[index]*dVz_dz);                    

            Tyy[index] += dt*(C12[index]*dVx_dx + C26[index]*dVx_dy + C25[index]*dVx_dz +
                              C26[index]*dVy_dx + C22[index]*dVy_dy + C24[index]*dVy_dz +
                              C25[index]*dVz_dx + C24[index]*dVz_dy + C23[index]*dVz_dz);                    

            Tzz[index] += dt*(C13[index]*dVx_dx + C36[index]*dVx_dy + C35[index]*dVx_dz +
                              C36[index]*dVy_dx + C23[index]*dVy_dy + C34[index]*dVy_dz +
                              C35[index]*dVz_dx + C34[index]*dVz_dy + C33[index]*dVz_dz);                                
        }

        if((i >= 3) && (i < nzz-4) && (j > 3) && (j < nxx-3) && (k > 3) && (k < nyy-3)) 
        {
            float dVx_dx = (FDM1*(Vx[i + (j-4)*nzz + k*nxx*nzz] - Vx[i + (j+3)*nzz + k*nxx*nzz]) +
                            FDM2*(Vx[i + (j+2)*nzz + k*nxx*nzz] - Vx[i + (j-3)*nzz + k*nxx*nzz]) +
                            FDM3*(Vx[i + (j-2)*nzz + k*nxx*nzz] - Vx[i + (j+1)*nzz + k*nxx*nzz]) +
                            FDM4*(Vx[i + j*nzz + k*nxx*nzz]     - Vx[i + (j-1)*nzz + k*nxx*nzz])) / dx;

            float dVy_dx = (FDM1*(Vy[i + (j-4)*nzz + k*nxx*nzz] - Vy[i + (j+3)*nzz + k*nxx*nzz]) +
                            FDM2*(Vy[i + (j+2)*nzz + k*nxx*nzz] - Vy[i + (j-3)*nzz + k*nxx*nzz]) +
                            FDM3*(Vy[i + (j-2)*nzz + k*nxx*nzz] - Vy[i + (j+1)*nzz + k*nxx*nzz]) +
                            FDM4*(Vy[i + j*nzz + k*nxx*nzz]     - Vy[i + (j-1)*nzz + k*nxx*nzz])) / dx;

            float dVz_dx = (FDM1*(Vz[i + (j-4)*nzz + k*nxx*nzz] - Vz[i + (j+3)*nzz + k*nxx*nzz]) +
                            FDM2*(Vz[i + (j+2)*nzz + k*nxx*nzz] - Vz[i + (j-3)*nzz + k*nxx*nzz]) +
                            FDM3*(Vz[i + (j-2)*nzz + k*nxx*nzz] - Vz[i + (j+1)*nzz + k*nxx*nzz]) +
                            FDM4*(Vz[i + j*nzz + k*nxx*nzz] -     Vz[i + (j-1)*nzz + k*nxx*nzz])) / dx;

            float dVx_dy = (FDM1*(Vx[i + j*nzz + (k-4)*nxx*nzz] - Vx[i + j*nzz + (k+3)*nxx*nzz]) +
                            FDM2*(Vx[i + j*nzz + (k+2)*nxx*nzz] - Vx[i + j*nzz + (k-3)*nxx*nzz]) +
                            FDM3*(Vx[i + j*nzz + (k-2)*nxx*nzz] - Vx[i + j*nzz + (k+1)*nxx*nzz]) +
                            FDM4*(Vx[i + j*nzz + k*nxx*nzz]     - Vx[i + j*nzz + (k-1)*nxx*nzz])) / dy;

            float dVy_dy = (FDM1*(Vy[i + j*nzz + (k-4)*nxx*nzz] - Vy[i + j*nzz + (k+3)*nxx*nzz]) +
                            FDM2*(Vy[i + j*nzz + (k+2)*nxx*nzz] - Vy[i + j*nzz + (k-3)*nxx*nzz]) +
                            FDM3*(Vy[i + j*nzz + (k-2)*nxx*nzz] - Vy[i + j*nzz + (k+1)*nxx*nzz]) +
                            FDM4*(Vy[i + j*nzz + k*nxx*nzz]     - Vy[i + j*nzz + (k-1)*nxx*nzz])) / dy;

            float dVz_dy = (FDM1*(Vz[i + j*nzz + (k-4)*nxx*nzz] - Vz[i + j*nzz + (k+3)*nxx*nzz]) +
                            FDM2*(Vz[i + j*nzz + (k+2)*nxx*nzz] - Vz[i + j*nzz + (k-3)*nxx*nzz]) +
                            FDM3*(Vz[i + j*nzz + (k-2)*nxx*nzz] - Vz[i + j*nzz + (k+1)*nxx*nzz]) +
                            FDM4*(Vz[i + j*nzz + k*nxx*nzz]     - Vz[i + j*nzz + (k-1)*nxx*nzz])) / dy;

            float dVx_dz = (FDM1*(Vx[(i-3) + j*nzz + k*nxx*nzz] - Vx[(i+4) + j*nzz + k*nxx*nzz]) +
                            FDM2*(Vx[(i+3) + j*nzz + k*nxx*nzz] - Vx[(i-2) + j*nzz + k*nxx*nzz]) +
                            FDM3*(Vx[(i-1) + j*nzz + k*nxx*nzz] - Vx[(i+2) + j*nzz + k*nxx*nzz]) +
                            FDM4*(Vx[(i+1) + j*nzz + k*nxx*nzz] - Vx[i + j*nzz + k*nxx*nzz])) / dz;

            float dVy_dz = (FDM1*(Vy[(i-3) + j*nzz + k*nxx*nzz] - Vy[(i+4) + j*nzz + k*nxx*nzz]) +
                            FDM2*(Vy[(i+3) + j*nzz + k*nxx*nzz] - Vy[(i-2) + j*nzz + k*nxx*nzz]) +
                            FDM3*(Vy[(i-1) + j*nzz + k*nxx*nzz] - Vy[(i+2) + j*nzz + k*nxx*nzz]) +
                            FDM4*(Vy[(i+1) + j*nzz + k*nxx*nzz] - Vy[i + j*nzz + k*nxx*nzz])) / dz;

            float dVz_dz = (FDM1*(Vz[(i-3) + j*nzz + k*nxx*nzz] - Vz[(i+4) + j*nzz + k*nxx*nzz]) +
                            FDM2*(Vz[(i+3) + j*nzz + k*nxx*nzz] - Vz[(i-2) + j*nzz + k*nxx*nzz]) +
                            FDM3*(Vz[(i-1) + j*nzz + k*nxx*nzz] - Vz[(i+2) + j*nzz + k*nxx*nzz]) +
                            FDM4*(Vz[(i+1) + j*nzz + k*nxx*nzz] - Vz[i + j*nzz + k*nxx*nzz])) / dz;
            
            float C16_xy = powf(0.25f*(1.0f/C16[i + (j+1)*nzz + (k+1)*nxx*nzz] + 1.0f/C16[i + (j+1)*nzz + k*nxx*nzz] + 
                                       1.0f/C16[i + j*nzz + (k+1)*nxx*nzz]     + 1.0f/C16[i + j*nzz + k*nxx*nzz]), -1.0f);

            float C26_xy = powf(0.25f*(1.0f/C26[i + (j+1)*nzz + (k+1)*nxx*nzz] + 1.0f/C26[i + (j+1)*nzz + k*nxx*nzz] + 
                                       1.0f/C26[i + j*nzz + (k+1)*nxx*nzz]     + 1.0f/C26[i + j*nzz + k*nxx*nzz]), -1.0f);

            float C36_xy = powf(0.25f*(1.0f/C36[i + (j+1)*nzz + (k+1)*nxx*nzz] + 1.0f/C36[i + (j+1)*nzz + k*nxx*nzz] + 
                                       1.0f/C36[i + j*nzz + (k+1)*nxx*nzz]     + 1.0f/C36[i + j*nzz + k*nxx*nzz]), -1.0f);

            float C46_xy = powf(0.25f*(1.0f/C46[i + (j+1)*nzz + (k+1)*nxx*nzz] + 1.0f/C46[i + (j+1)*nzz + k*nxx*nzz] + 
                                       1.0f/C46[i + j*nzz + (k+1)*nxx*nzz]     + 1.0f/C46[i + j*nzz + k*nxx*nzz]), -1.0f);

            float C56_xy = powf(0.25f*(1.0f/C56[i + (j+1)*nzz + (k+1)*nxx*nzz] + 1.0f/C56[i + (j+1)*nzz + k*nxx*nzz] + 
                                       1.0f/C56[i + j*nzz + (k+1)*nxx*nzz]     + 1.0f/C56[i + j*nzz + k*nxx*nzz]), -1.0f);

            float C66_xy = powf(0.25f*(1.0f/C66[i + (j+1)*nzz + (k+1)*nxx*nzz] + 1.0f/C66[i + (j+1)*nzz + k*nxx*nzz] + 
                                       1.0f/C66[i + j*nzz + (k+1)*nxx*nzz]     + 1.0f/C66[i + j*nzz + k*nxx*nzz]), -1.0f);

            Txy[index] += dt*(C16_xy*dVx_dx + C66_xy*dVx_dy + C56_xy*dVx_dz +
                              C66_xy*dVy_dx + C26_xy*dVy_dy + C46_xy*dVy_dz +
                              C56_xy*dVz_dx + C46_xy*dVz_dy + C36_xy*dVz_dz);                    
        }

        if((i > 3) && (i < nzz-3) && (j > 3) && (j < nxx-3) && (k >= 3) && (k < nyy-4)) 
        {
            float dVx_dx = (FDM1*(Vx[i + (j-4)*nzz + k*nxx*nzz] - Vx[i + (j+3)*nzz + k*nxx*nzz]) +
                            FDM2*(Vx[i + (j+2)*nzz + k*nxx*nzz] - Vx[i + (j-3)*nzz + k*nxx*nzz]) +
                            FDM3*(Vx[i + (j-2)*nzz + k*nxx*nzz] - Vx[i + (j+1)*nzz + k*nxx*nzz]) +
                            FDM4*(Vx[i + j*nzz + k*nxx*nzz]     - Vx[i + (j-1)*nzz + k*nxx*nzz])) / dx;

            float dVy_dx = (FDM1*(Vy[i + (j-4)*nzz + k*nxx*nzz] - Vy[i + (j+3)*nzz + k*nxx*nzz]) +
                            FDM2*(Vy[i + (j+2)*nzz + k*nxx*nzz] - Vy[i + (j-3)*nzz + k*nxx*nzz]) +
                            FDM3*(Vy[i + (j-2)*nzz + k*nxx*nzz] - Vy[i + (j+1)*nzz + k*nxx*nzz]) +
                            FDM4*(Vy[i + j*nzz + k*nxx*nzz]     - Vy[i + (j-1)*nzz + k*nxx*nzz])) / dx;

            float dVz_dx = (FDM1*(Vz[i + (j-4)*nzz + k*nxx*nzz] - Vz[i + (j+3)*nzz + k*nxx*nzz]) +
                            FDM2*(Vz[i + (j+2)*nzz + k*nxx*nzz] - Vz[i + (j-3)*nzz + k*nxx*nzz]) +
                            FDM3*(Vz[i + (j-2)*nzz + k*nxx*nzz] - Vz[i + (j+1)*nzz + k*nxx*nzz]) +
                            FDM4*(Vz[i + j*nzz + k*nxx*nzz]     - Vz[i + (j-1)*nzz + k*nxx*nzz])) / dx;

            float dVx_dy = (FDM1*(Vx[i + j*nzz + (k-3)*nxx*nzz] - Vx[i + j*nzz + (k+4)*nxx*nzz]) +
                            FDM2*(Vx[i + j*nzz + (k+3)*nxx*nzz] - Vx[i + j*nzz + (k-2)*nxx*nzz]) +
                            FDM3*(Vx[i + j*nzz + (k-1)*nxx*nzz] - Vx[i + j*nzz + (k+2)*nxx*nzz]) +
                            FDM4*(Vx[i + j*nzz + (k+1)*nxx*nzz] - Vx[i + j*nzz + k*nxx*nzz])) / dy;

            float dVy_dy = (FDM1*(Vy[i + j*nzz + (k-3)*nxx*nzz] - Vy[i + j*nzz + (k+4)*nxx*nzz]) +
                            FDM2*(Vy[i + j*nzz + (k+3)*nxx*nzz] - Vy[i + j*nzz + (k-2)*nxx*nzz]) +
                            FDM3*(Vy[i + j*nzz + (k-1)*nxx*nzz] - Vy[i + j*nzz + (k+2)*nxx*nzz]) +
                            FDM4*(Vy[i + j*nzz + (k+1)*nxx*nzz] - Vy[i + j*nzz + k*nxx*nzz])) / dy;

            float dVz_dy = (FDM1*(Vz[i + j*nzz + (k-3)*nxx*nzz] - Vz[i + j*nzz + (k+4)*nxx*nzz]) +
                            FDM2*(Vz[i + j*nzz + (k+3)*nxx*nzz] - Vz[i + j*nzz + (k-2)*nxx*nzz]) +
                            FDM3*(Vz[i + j*nzz + (k-1)*nxx*nzz] - Vz[i + j*nzz + (k+2)*nxx*nzz]) +
                            FDM4*(Vz[i + j*nzz + (k+1)*nxx*nzz] - Vz[i + j*nzz + k*nxx*nzz])) / dy;

            float dVx_dz = (FDM1*(Vx[(i-4) + j*nzz + k*nxx*nzz] - Vx[(i+3) + j*nzz + k*nxx*nzz]) +
                            FDM2*(Vx[(i+2) + j*nzz + k*nxx*nzz] - Vx[(i-3) + j*nzz + k*nxx*nzz]) +
                            FDM3*(Vx[(i-2) + j*nzz + k*nxx*nzz] - Vx[(i+1) + j*nzz + k*nxx*nzz]) +
                            FDM4*(Vx[i + j*nzz + k*nxx*nzz]     - Vx[(i-1) + j*nzz + k*nxx*nzz])) / dz;

            float dVy_dz = (FDM1*(Vy[(i+4) + j*nzz + k*nxx*nzz] - Vy[(i+3) + j*nzz + k*nxx*nzz]) +
                            FDM2*(Vy[(i+2) + j*nzz + k*nxx*nzz] - Vy[(i-3) + j*nzz + k*nxx*nzz]) +
                            FDM3*(Vy[(i-2) + j*nzz + k*nxx*nzz] - Vy[(i+1) + j*nzz + k*nxx*nzz]) +
                            FDM4*(Vy[i + j*nzz + k*nxx*nzz]     - Vy[(i-1) + j*nzz + k*nxx*nzz])) / dz;

            float dVz_dz = (FDM1*(Vz[(i-4) + j*nzz + k*nxx*nzz] - Vz[(i+3) + j*nzz + k*nxx*nzz]) +
                            FDM2*(Vz[(i+2) + j*nzz + k*nxx*nzz] - Vz[(i-3) + j*nzz + k*nxx*nzz]) +
                            FDM3*(Vz[(i-2) + j*nzz + k*nxx*nzz] - Vz[(i+1) + j*nzz + k*nxx*nzz]) +
                            FDM4*(Vz[i + j*nzz + k*nxx*nzz]     - Vz[(i-1) + j*nzz + k*nxx*nzz])) / dz;

            float C15_xz = powf(0.25f*(1.0f/C15[(i+1) + (j+1)*nzz + k*nxx*nzz] + 1.0f/C15[i + (j+1)*nzz + k*nxx*nzz] + 
                                       1.0f/C15[(i+1) + j*nzz + k*nxx*nzz]     + 1.0f/C15[i + j*nzz + k*nxx*nzz]), -1.0f);

            float C25_xz = powf(0.25f*(1.0f/C25[(i+1) + (j+1)*nzz + k*nxx*nzz] + 1.0f/C25[i + (j+1)*nzz + k*nxx*nzz] + 
                                       1.0f/C25[(i+1) + j*nzz + k*nxx*nzz]     + 1.0f/C25[i + j*nzz + k*nxx*nzz]), -1.0f);

            float C35_xz = powf(0.25f*(1.0f/C35[(i+1) + (j+1)*nzz + k*nxx*nzz] + 1.0f/C35[i + (j+1)*nzz + k*nxx*nzz] + 
                                       1.0f/C35[(i+1) + j*nzz + k*nxx*nzz]     + 1.0f/C35[i + j*nzz + k*nxx*nzz]), -1.0f);

            float C45_xz = powf(0.25f*(1.0f/C45[(i+1) + (j+1)*nzz + k*nxx*nzz] + 1.0f/C45[i + (j+1)*nzz + k*nxx*nzz] + 
                                       1.0f/C45[(i+1) + j*nzz + k*nxx*nzz]     + 1.0f/C45[i + j*nzz + k*nxx*nzz]), -1.0f);

            float C55_xz = powf(0.25f*(1.0f/C55[(i+1) + (j+1)*nzz + k*nxx*nzz] + 1.0f/C55[i + (j+1)*nzz + k*nxx*nzz] + 
                                       1.0f/C55[(i+1) + j*nzz + k*nxx*nzz]     + 1.0f/C55[i + j*nzz + k*nxx*nzz]), -1.0f);

            float C56_xz = powf(0.25f*(1.0f/C56[(i+1) + (j+1)*nzz + k*nxx*nzz] + 1.0f/C56[i + (j+1)*nzz + k*nxx*nzz] + 
                                       1.0f/C56[(i+1) + j*nzz + k*nxx*nzz]     + 1.0f/C56[i + j*nzz + k*nxx*nzz]), -1.0f);

            Txz[index] += dt*(C15_xz*dVx_dx + C56_xz*dVx_dy + C55_xz*dVx_dz +
                              C56_xz*dVy_dx + C25_xz*dVy_dy + C45_xz*dVy_dz +
                              C55_xz*dVz_dx + C45_xz*dVz_dy + C35_xz*dVz_dz);                    
        }

        if((i > 3) && (i < nzz-3) && (j >= 3) && (j < nxx-4) && (k > 3) && (k < nyy-3)) 
        {
            float dVx_dx = (FDM1*(Vx[i + (j-3)*nzz + k*nxx*nzz] - Vx[i + (j+4)*nzz + k*nxx*nzz]) +
                            FDM2*(Vx[i + (j+3)*nzz + k*nxx*nzz] - Vx[i + (j-2)*nzz + k*nxx*nzz]) +
                            FDM3*(Vx[i + (j-1)*nzz + k*nxx*nzz] - Vx[i + (j+2)*nzz + k*nxx*nzz]) +
                            FDM4*(Vx[i + (j+1)*nzz + k*nxx*nzz] - Vx[i + j*nzz + k*nxx*nzz])) / dx;

            float dVy_dx = (FDM1*(Vy[i + (j-3)*nzz + k*nxx*nzz] - Vy[i + (j+4)*nzz + k*nxx*nzz]) +
                            FDM2*(Vy[i + (j+3)*nzz + k*nxx*nzz] - Vy[i + (j-2)*nzz + k*nxx*nzz]) +
                            FDM3*(Vy[i + (j-1)*nzz + k*nxx*nzz] - Vy[i + (j+2)*nzz + k*nxx*nzz]) +
                            FDM4*(Vy[i + (j+1)*nzz + k*nxx*nzz] - Vy[i + j*nzz + k*nxx*nzz])) / dx;

            float dVz_dx = (FDM1*(Vz[i + (j-3)*nzz + k*nxx*nzz] - Vz[i + (j+4)*nzz + k*nxx*nzz]) +
                            FDM2*(Vz[i + (j+3)*nzz + k*nxx*nzz] - Vz[i + (j-2)*nzz + k*nxx*nzz]) +
                            FDM3*(Vz[i + (j-1)*nzz + k*nxx*nzz] - Vz[i + (j+2)*nzz + k*nxx*nzz]) +
                            FDM4*(Vz[i + (j+1)*nzz + k*nxx*nzz] - Vz[i + j*nzz + k*nxx*nzz])) / dx;

            float dVx_dy = (FDM1*(Vx[i + j*nzz + (k-4)*nxx*nzz] - Vx[i + j*nzz + (k+3)*nxx*nzz]) +
                            FDM2*(Vx[i + j*nzz + (k+2)*nxx*nzz] - Vx[i + j*nzz + (k-3)*nxx*nzz]) +
                            FDM3*(Vx[i + j*nzz + (k-2)*nxx*nzz] - Vx[i + j*nzz + (k+1)*nxx*nzz]) +
                            FDM4*(Vx[i + j*nzz + k*nxx*nzz]     - Vx[i + j*nzz + (k-1)*nxx*nzz])) / dy;

            float dVy_dy = (FDM1*(Vy[i + j*nzz + (k-4)*nxx*nzz] - Vy[i + j*nzz + (k+3)*nxx*nzz]) +
                            FDM2*(Vy[i + j*nzz + (k+2)*nxx*nzz] - Vy[i + j*nzz + (k-3)*nxx*nzz]) +
                            FDM3*(Vy[i + j*nzz + (k-2)*nxx*nzz] - Vy[i + j*nzz + (k+1)*nxx*nzz]) +
                            FDM4*(Vy[i + j*nzz + k*nxx*nzz]     - Vy[i + j*nzz + (k-1)*nxx*nzz])) / dy;

            float dVz_dy = (FDM1*(Vz[i + j*nzz + (k-4)*nxx*nzz] - Vz[i + j*nzz + (k+3)*nxx*nzz]) +
                            FDM2*(Vz[i + j*nzz + (k+2)*nxx*nzz] - Vz[i + j*nzz + (k-3)*nxx*nzz]) +
                            FDM3*(Vz[i + j*nzz + (k-2)*nxx*nzz] - Vz[i + j*nzz + (k+1)*nxx*nzz]) +
                            FDM4*(Vz[i + j*nzz + k*nxx*nzz]     - Vz[i + j*nzz + (k-1)*nxx*nzz])) / dy;

            float dVx_dz = (FDM1*(Vx[(i-4) + j*nzz + k*nxx*nzz] - Vx[(i+3) + j*nzz + k*nxx*nzz]) +
                            FDM2*(Vx[(i+2) + j*nzz + k*nxx*nzz] - Vx[(i-3) + j*nzz + k*nxx*nzz]) +
                            FDM3*(Vx[(i-2) + j*nzz + k*nxx*nzz] - Vx[(i+1) + j*nzz + k*nxx*nzz]) +
                            FDM4*(Vx[i + j*nzz + k*nxx*nzz]     - Vx[(i-1) + j*nzz + k*nxx*nzz])) / dz;

            float dVy_dz = (FDM1*(Vy[(i-4) + j*nzz + k*nxx*nzz] - Vy[(i+3) + j*nzz + k*nxx*nzz]) +
                            FDM2*(Vy[(i+2) + j*nzz + k*nxx*nzz] - Vy[(i-3) + j*nzz + k*nxx*nzz]) +
                            FDM3*(Vy[(i-2) + j*nzz + k*nxx*nzz] - Vy[(i+1) + j*nzz + k*nxx*nzz]) +
                            FDM4*(Vy[i + j*nzz + k*nxx*nzz]     - Vy[(i-1) + j*nzz + k*nxx*nzz])) / dz;

            float dVz_dz = (FDM1*(Vz[(i-4) + j*nzz + k*nxx*nzz] - Vz[(i+3) + j*nzz + k*nxx*nzz]) +
                            FDM2*(Vz[(i+2) + j*nzz + k*nxx*nzz] - Vz[(i-3) + j*nzz + k*nxx*nzz]) +
                            FDM3*(Vz[(i-2) + j*nzz + k*nxx*nzz] - Vz[(i+1) + j*nzz + k*nxx*nzz]) +
                            FDM4*(Vz[i + j*nzz + k*nxx*nzz]     - Vz[(i-1) + j*nzz + k*nxx*nzz])) / dz;

            float C14_yz = powf(0.25f*(1.0f/C14[(i+1) + j*nzz + (k+1)*nxx*nzz] + 1.0f/C14[i + j*nzz + (k+1)*nxx*nzz] + 
                                       1.0f/C14[(i+1) + j*nzz + k*nxx*nzz] +     1.0f/C14[i + j*nzz + k*nxx*nzz]), -1.0f);

            float C24_yz = powf(0.25f*(1.0f/C24[(i+1) + j*nzz + (k+1)*nxx*nzz] + 1.0f/C24[i + j*nzz + (k+1)*nxx*nzz] + 
                                       1.0f/C24[(i+1) + j*nzz + k*nxx*nzz] +     1.0f/C24[i + j*nzz + k*nxx*nzz]), -1.0f);

            float C34_yz = powf(0.25f*(1.0f/C34[(i+1) + j*nzz + (k+1)*nxx*nzz] + 1.0f/C34[i + j*nzz + (k+1)*nxx*nzz] + 
                                       1.0f/C34[(i+1) + j*nzz + k*nxx*nzz] +     1.0f/C34[i + j*nzz + k*nxx*nzz]), -1.0f);

            float C44_yz = powf(0.25f*(1.0f/C44[(i+1) + j*nzz + (k+1)*nxx*nzz] + 1.0f/C44[i + j*nzz + (k+1)*nxx*nzz] + 
                                       1.0f/C44[(i+1) + j*nzz + k*nxx*nzz] +     1.0f/C44[i + j*nzz + k*nxx*nzz]), -1.0f);

            float C45_yz = powf(0.25f*(1.0f/C45[(i+1) + j*nzz + (k+1)*nxx*nzz] + 1.0f/C45[i + j*nzz + (k+1)*nxx*nzz] + 
                                       1.0f/C45[(i+1) + j*nzz + k*nxx*nzz] +     1.0f/C45[i + j*nzz + k*nxx*nzz]), -1.0f);

            float C46_yz = powf(0.25f*(1.0f/C46[(i+1) + j*nzz + (k+1)*nxx*nzz] + 1.0f/C46[i + j*nzz + (k+1)*nxx*nzz] + 
                                       1.0f/C46[(i+1) + j*nzz + k*nxx*nzz] +     1.0f/C46[i + j*nzz + k*nxx*nzz]), -1.0f);

            Tyz[index] += dt*(C14_yz*dVx_dx + C46_yz*dVx_dy + C45_yz*dVx_dz +
                              C46_yz*dVy_dx + C24_yz*dVy_dy + C44_yz*dVy_dz +
                              C45_yz*dVz_dx + C44_yz*dVz_dy + C34_yz*dVz_dz);                    
        }

        if ((i > 3) && (i < nzz-4) && (j > 3) && (j < nxx-4) && (k > 3) && (k < nyy-4))
        {
            P[index] = (Txx[index] + Tyy[index] + Tzz[index]) / 3.0f;
        }
    }
}

