# include "elastic_ani.cuh"

void Elastic_ANI::set_conditions()
{
    modeling_type = "elastic_ani";
    modeling_name = "Modeling type: Elastic anisotropic solver";

    eikonal = new Eikonal_ANI();
    eikonal->parameters = parameters;
    eikonal->set_parameters();

    int samples = 2*RSGR+1;
    int nKernel = samples*samples*samples;
    float * kernel = new float[nKernel]();

    int index = 0;
    float sum = 0.0f;

    for (int z = -RSGR; z <= RSGR; z++)
    {
        for (int x = -RSGR; x <= RSGR; x++)
        {
            for (int y = -RSGR; y <= RSGR; y++)
            {          
                float r = sqrtf(x*x + y*y + z*z);

                kernel[index] = 1.0f/sqrtf(2.0f*M_PI)*expf(-0.5f*r*r);
    
                sum += kernel[index]; 
            
                ++index;
            }
        }
    }

    for (index = 0; index < nKernel; index++) 
        kernel[index] /= sum;

    cudaMalloc((void**)&(dwc), nKernel*sizeof(float));
    cudaMemcpy(dwc, kernel, nKernel*sizeof(float), cudaMemcpyHostToDevice);
    
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
        compute_velocity_rsg<<<nBlocks, nThreads>>>(d_Vx, d_Vy, d_Vz, d_Txx, d_Tyy, d_Tzz, d_Txz, d_Tyz, d_Txy, d_B, d_T, d1D, d2D, d3D, wavelet, dwc, dx, dy, dz, dt, tId, tlag, sIdx, sIdy, sIdz, nxx, nyy, nzz, nb, nt);
        cudaDeviceSynchronize();

        compute_pressure_rsg<<<nBlocks, nThreads>>>(d_Vx, d_Vy, d_Vz, d_Txx, d_Tyy, d_Tzz, d_Txz, d_Tyz, d_Txy, d_P, d_T, d_C11, d_C12, d_C13, d_C14, d_C15, d_C16, d_C22, d_C23, d_C24, d_C25, d_C26, d_C33, d_C34, d_C35, d_C36, d_C44, d_C45, d_C46, d_C55, d_C56, d_C66, tId, tlag, dx, dy, dz, dt, nxx, nyy, nzz);
        cudaDeviceSynchronize();

        compute_seismogram<<<sBlocks, nThreads>>>(d_P, rIdx, rIdy, rIdz, seismogram, geometry->spread[srcId], tId, tlag, nt, nxx, nzz);     
        cudaDeviceSynchronize();
    }
}

__global__ void compute_velocity_rsg(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * B, float * T,  
                                     float * damp1D, float * damp2D, float * damp3D, float * wavelet, float * dwc, float dx, float dy, float dz, float dt, int tId, int tlag, 
                                     int sIdx, int sIdy, int sIdz, int nxx, int nyy, int nzz, int nb, int nt)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         
    int j = (int) (index - k*nxx*nzz) / nzz;   
    int i = (int) (index - j*nzz - k*nxx*nzz); 

    if ((index == 0) && (tId < nt))
    {
        int sId = 0;

        for (int si = -RSGR; si <= RSGR; si++)
        {
            for (int sj = -RSGR; sj <= RSGR; sj++)
            {
                for (int sk = -RSGR; sk <= RSGR; sk++)
                {
                    Txx[(sIdz+si) + (sIdx+sj)*nzz + (sIdy+sk)*nxx*nzz] += dwc[sId]*wavelet[tId] / (dx*dy*dz);
                    Tyy[(sIdz+si) + (sIdx+sj)*nzz + (sIdy+sk)*nxx*nzz] += dwc[sId]*wavelet[tId] / (dx*dy*dz);
                    Tzz[(sIdz+si) + (sIdx+sj)*nzz + (sIdy+sk)*nxx*nzz] += dwc[sId]*wavelet[tId] / (dx*dy*dz);                            
                    
                    ++sId;
                }
            }
        }
    }

    float d1_Txx = 0.0f; float d2_Txx = 0.0f; float d3_Txx = 0.0f; float d4_Txx = 0.0f;
    float d1_Tyy = 0.0f; float d2_Tyy = 0.0f; float d3_Tyy = 0.0f; float d4_Tyy = 0.0f;
    float d1_Tzz = 0.0f; float d2_Tzz = 0.0f; float d3_Tzz = 0.0f; float d4_Tzz = 0.0f;
    float d1_Txy = 0.0f; float d2_Txy = 0.0f; float d3_Txy = 0.0f; float d4_Txy = 0.0f;
    float d1_Txz = 0.0f; float d2_Txz = 0.0f; float d3_Txz = 0.0f; float d4_Txz = 0.0f;
    float d1_Tyz = 0.0f; float d2_Tyz = 0.0f; float d3_Tyz = 0.0f; float d4_Tyz = 0.0f;            
 
    float FDM[] = {FDM4, -FDM3, FDM2, -FDM1};
    
    if ((T[index] < (float)(tId + tlag)*dt) && (index < nxx*nyy*nzz))
    {
        if((i >= 3) && (i < nzz-4) && (j >= 3) && (j < nxx-4) && (k >= 3) && (k < nyy-4)) 
        {    
            for (int rsg = 0; rsg < 4; rsg++)
            {
                d1_Txx += FDM[rsg]*(Txx[(i+rsg+1) + (j+rsg+1)*nzz + (k+rsg+1)*nxx*nzz] - Txx[(i-rsg) + (j-rsg)*nzz + (k-rsg)*nxx*nzz]);
                d1_Tyy += FDM[rsg]*(Tyy[(i+rsg+1) + (j+rsg+1)*nzz + (k+rsg+1)*nxx*nzz] - Tyy[(i-rsg) + (j-rsg)*nzz + (k-rsg)*nxx*nzz]);
                d1_Tzz += FDM[rsg]*(Tzz[(i+rsg+1) + (j+rsg+1)*nzz + (k+rsg+1)*nxx*nzz] - Tzz[(i-rsg) + (j-rsg)*nzz + (k-rsg)*nxx*nzz]);
                d1_Txy += FDM[rsg]*(Txy[(i+rsg+1) + (j+rsg+1)*nzz + (k+rsg+1)*nxx*nzz] - Txy[(i-rsg) + (j-rsg)*nzz + (k-rsg)*nxx*nzz]);
                d1_Txz += FDM[rsg]*(Txz[(i+rsg+1) + (j+rsg+1)*nzz + (k+rsg+1)*nxx*nzz] - Txz[(i-rsg) + (j-rsg)*nzz + (k-rsg)*nxx*nzz]);
                d1_Tyz += FDM[rsg]*(Tyz[(i+rsg+1) + (j+rsg+1)*nzz + (k+rsg+1)*nxx*nzz] - Tyz[(i-rsg) + (j-rsg)*nzz + (k-rsg)*nxx*nzz]);

                d2_Txx += FDM[rsg]*(Txx[(i-rsg) + (j+rsg+1)*nzz + (k+rsg+1)*nxx*nzz] - Txx[(i+rsg+1) + (j-rsg)*nzz + (k-rsg)*nxx*nzz]);
                d2_Tyy += FDM[rsg]*(Tyy[(i-rsg) + (j+rsg+1)*nzz + (k+rsg+1)*nxx*nzz] - Tyy[(i+rsg+1) + (j-rsg)*nzz + (k-rsg)*nxx*nzz]);
                d2_Tzz += FDM[rsg]*(Tzz[(i-rsg) + (j+rsg+1)*nzz + (k+rsg+1)*nxx*nzz] - Tzz[(i+rsg+1) + (j-rsg)*nzz + (k-rsg)*nxx*nzz]);
                d2_Txy += FDM[rsg]*(Txy[(i-rsg) + (j+rsg+1)*nzz + (k+rsg+1)*nxx*nzz] - Txy[(i+rsg+1) + (j-rsg)*nzz + (k-rsg)*nxx*nzz]);
                d2_Txz += FDM[rsg]*(Txz[(i-rsg) + (j+rsg+1)*nzz + (k+rsg+1)*nxx*nzz] - Txz[(i+rsg+1) + (j-rsg)*nzz + (k-rsg)*nxx*nzz]);
                d2_Tyz += FDM[rsg]*(Tyz[(i-rsg) + (j+rsg+1)*nzz + (k+rsg+1)*nxx*nzz] - Tyz[(i+rsg+1) + (j-rsg)*nzz + (k-rsg)*nxx*nzz]);

                d3_Txx += FDM[rsg]*(Txx[(i+rsg+1) + (j+rsg+1)*nzz + (k-rsg)*nxx*nzz] - Txx[(i-rsg) + (j-rsg)*nzz + (k+rsg+1)*nxx*nzz]);
                d3_Tyy += FDM[rsg]*(Tyy[(i+rsg+1) + (j+rsg+1)*nzz + (k-rsg)*nxx*nzz] - Tyy[(i-rsg) + (j-rsg)*nzz + (k+rsg+1)*nxx*nzz]);
                d3_Tzz += FDM[rsg]*(Tzz[(i+rsg+1) + (j+rsg+1)*nzz + (k-rsg)*nxx*nzz] - Tzz[(i-rsg) + (j-rsg)*nzz + (k+rsg+1)*nxx*nzz]);
                d3_Txy += FDM[rsg]*(Txy[(i+rsg+1) + (j+rsg+1)*nzz + (k-rsg)*nxx*nzz] - Txy[(i-rsg) + (j-rsg)*nzz + (k+rsg+1)*nxx*nzz]);
                d3_Txz += FDM[rsg]*(Txz[(i+rsg+1) + (j+rsg+1)*nzz + (k-rsg)*nxx*nzz] - Txz[(i-rsg) + (j-rsg)*nzz + (k+rsg+1)*nxx*nzz]);
                d3_Tyz += FDM[rsg]*(Tyz[(i+rsg+1) + (j+rsg+1)*nzz + (k-rsg)*nxx*nzz] - Tyz[(i-rsg) + (j-rsg)*nzz + (k+rsg+1)*nxx*nzz]);

                d4_Txx += FDM[rsg]*(Txx[(i-rsg) + (j+rsg+1)*nzz + (k-rsg)*nxx*nzz] - Txx[(i+rsg+1) + (j-rsg)*nzz + (k+rsg+1)*nxx*nzz]);
                d4_Tyy += FDM[rsg]*(Tyy[(i-rsg) + (j+rsg+1)*nzz + (k-rsg)*nxx*nzz] - Tyy[(i+rsg+1) + (j-rsg)*nzz + (k+rsg+1)*nxx*nzz]);
                d4_Tzz += FDM[rsg]*(Tzz[(i-rsg) + (j+rsg+1)*nzz + (k-rsg)*nxx*nzz] - Tzz[(i+rsg+1) + (j-rsg)*nzz + (k+rsg+1)*nxx*nzz]);
                d4_Txy += FDM[rsg]*(Txy[(i-rsg) + (j+rsg+1)*nzz + (k-rsg)*nxx*nzz] - Txy[(i+rsg+1) + (j-rsg)*nzz + (k+rsg+1)*nxx*nzz]);
                d4_Txz += FDM[rsg]*(Txz[(i-rsg) + (j+rsg+1)*nzz + (k-rsg)*nxx*nzz] - Txz[(i+rsg+1) + (j-rsg)*nzz + (k+rsg+1)*nxx*nzz]);
                d4_Tyz += FDM[rsg]*(Tyz[(i-rsg) + (j+rsg+1)*nzz + (k-rsg)*nxx*nzz] - Tyz[(i+rsg+1) + (j-rsg)*nzz + (k+rsg+1)*nxx*nzz]);
            }
        }
    
        float dTxx_dx = 0.25f*(d1_Txx + d2_Txx + d3_Txx + d4_Txx) / dx;
        float dTxy_dx = 0.25f*(d1_Txy + d2_Txy + d3_Txy + d4_Txy) / dx;
        float dTxz_dx = 0.25f*(d1_Txz + d2_Txz + d3_Txz + d4_Txz) / dx;
    
        float dTxy_dy = 0.25f*(d1_Txy + d2_Txy - d3_Txy - d4_Txy) / dy;
        float dTyy_dy = 0.25f*(d1_Tyy + d2_Tyy - d3_Tyy - d4_Tyy) / dy;
        float dTyz_dy = 0.25f*(d1_Tyz + d2_Tyz - d3_Tyz - d4_Tyz) / dy;
        
        float dTxz_dz = 0.25f*(d1_Txz - d2_Txz + d3_Txz - d4_Txz) / dz;
        float dTyz_dz = 0.25f*(d1_Tyz - d2_Tyz + d3_Tyz - d4_Tyz) / dz;
        float dTzz_dz = 0.25f*(d1_Tzz - d2_Tzz + d3_Tzz - d4_Tzz) / dz;
    
        Vx[index] += dt*B[index]*(dTxx_dx + dTxy_dy + dTxz_dz); 
        Vy[index] += dt*B[index]*(dTxy_dx + dTyy_dy + dTyz_dz);
        Vz[index] += dt*B[index]*(dTxz_dx + dTyz_dy + dTzz_dz);    
        
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

__global__ void compute_pressure_rsg(float * Vx, float * Vy, float * Vz, float * Txx, float * Tyy, float * Tzz, float * Txz, float * Tyz, float * Txy, float * P, float * T, 
                                     float * C11, float * C12, float * C13, float * C14, float * C15, float * C16, float * C22, float * C23, float * C24, float * C25, 
                                     float * C26, float * C33, float * C34, float * C35, float * C36, float * C44, float * C45, float * C46, float * C55, float * C56, 
                                     float * C66, int tId, int tlag, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         
    int j = (int) (index - k*nxx*nzz) / nzz;   
    int i = (int) (index - j*nzz - k*nxx*nzz); 

    float d1_Vx = 0.0f; float d2_Vx = 0.0f; float d3_Vx = 0.0f; float d4_Vx = 0.0f;
    float d1_Vy = 0.0f; float d2_Vy = 0.0f; float d3_Vy = 0.0f; float d4_Vy = 0.0f;
    float d1_Vz = 0.0f; float d2_Vz = 0.0f; float d3_Vz = 0.0f; float d4_Vz = 0.0f;

    float FDM[] = {FDM4, -FDM3, FDM2, -FDM1};

    if ((T[index] < (float)(tId + tlag)*dt) && (index < nxx*nyy*nzz))
    {
        if((i > 3) && (i < nzz-3) && (j > 3) && (j < nxx-3) && (k > 3) && (k < nyy-3)) 
        {
            for (int rsg = 0; rsg < 4; rsg++)
            {       
                d1_Vx += FDM[rsg]*(Vx[(i+rsg) + (j+rsg)*nzz + (k+rsg)*nxx*nzz] - Vx[(i-rsg-1) + (j-rsg-1)*nzz + (k-rsg-1)*nxx*nzz]);      
                d1_Vy += FDM[rsg]*(Vy[(i+rsg) + (j+rsg)*nzz + (k+rsg)*nxx*nzz] - Vy[(i-rsg-1) + (j-rsg-1)*nzz + (k-rsg-1)*nxx*nzz]);      
                d1_Vz += FDM[rsg]*(Vz[(i+rsg) + (j+rsg)*nzz + (k+rsg)*nxx*nzz] - Vz[(i-rsg-1) + (j-rsg-1)*nzz + (k-rsg-1)*nxx*nzz]);      
    
                d2_Vx += FDM[rsg]*(Vx[(i-rsg-1) + (j+rsg)*nzz + (k+rsg)*nxx*nzz] - Vx[(i+rsg) + (j-rsg-1)*nzz + (k-rsg-1)*nxx*nzz]);      
                d2_Vy += FDM[rsg]*(Vy[(i-rsg-1) + (j+rsg)*nzz + (k+rsg)*nxx*nzz] - Vy[(i+rsg) + (j-rsg-1)*nzz + (k-rsg-1)*nxx*nzz]);      
                d2_Vz += FDM[rsg]*(Vz[(i-rsg-1) + (j+rsg)*nzz + (k+rsg)*nxx*nzz] - Vz[(i+rsg) + (j-rsg-1)*nzz + (k-rsg-1)*nxx*nzz]);      
    
                d3_Vx += FDM[rsg]*(Vx[(i+rsg) + (j+rsg)*nzz + (k-rsg-1)*nxx*nzz] - Vx[(i-rsg-1) + (j-rsg-1)*nzz + (k+rsg)*nxx*nzz]);      
                d3_Vy += FDM[rsg]*(Vy[(i+rsg) + (j+rsg)*nzz + (k-rsg-1)*nxx*nzz] - Vy[(i-rsg-1) + (j-rsg-1)*nzz + (k+rsg)*nxx*nzz]);      
                d3_Vz += FDM[rsg]*(Vz[(i+rsg) + (j+rsg)*nzz + (k-rsg-1)*nxx*nzz] - Vz[(i-rsg-1) + (j-rsg-1)*nzz + (k+rsg)*nxx*nzz]);      
    
                d4_Vx += FDM[rsg]*(Vx[(i-rsg-1) + (j+rsg)*nzz + (k-rsg-1)*nxx*nzz] - Vx[(i+rsg) + (j-rsg-1)*nzz + (k+rsg)*nxx*nzz]);      
                d4_Vy += FDM[rsg]*(Vy[(i-rsg-1) + (j+rsg)*nzz + (k-rsg-1)*nxx*nzz] - Vy[(i+rsg) + (j-rsg-1)*nzz + (k+rsg)*nxx*nzz]);      
                d4_Vz += FDM[rsg]*(Vz[(i-rsg-1) + (j+rsg)*nzz + (k-rsg-1)*nxx*nzz] - Vz[(i+rsg) + (j-rsg-1)*nzz + (k+rsg)*nxx*nzz]);                          
            }
        }
    
        float dVx_dx = 0.25f*(d1_Vx + d2_Vx + d3_Vx + d4_Vx) / dx;
        float dVy_dx = 0.25f*(d1_Vy + d2_Vy + d3_Vy + d4_Vy) / dx;
        float dVz_dx = 0.25f*(d1_Vz + d2_Vz + d3_Vz + d4_Vz) / dx;
    
        float dVx_dy = 0.25f*(d1_Vx + d2_Vx - d3_Vx - d4_Vx) / dy;
        float dVy_dy = 0.25f*(d1_Vy + d2_Vy - d3_Vy - d4_Vy) / dy;
        float dVz_dy = 0.25f*(d1_Vz + d2_Vz - d3_Vz - d4_Vz) / dy;
    
        float dVx_dz = 0.25f*(d1_Vx - d2_Vx + d3_Vx - d4_Vx) / dz;
        float dVy_dz = 0.25f*(d1_Vy - d2_Vy + d3_Vy - d4_Vy) / dz;
        float dVz_dz = 0.25f*(d1_Vz - d2_Vz + d3_Vz - d4_Vz) / dz;
    
        Txx[index] += dt*(C11[index]*dVx_dx + C16[index]*dVx_dy + C15[index]*dVx_dz +
                          C16[index]*dVy_dx + C12[index]*dVy_dy + C14[index]*dVy_dz +
                          C15[index]*dVz_dx + C14[index]*dVz_dy + C13[index]*dVz_dz);                    
    
        Tyy[index] += dt*(C12[index]*dVx_dx + C26[index]*dVx_dy + C25[index]*dVx_dz +
                          C26[index]*dVy_dx + C22[index]*dVy_dy + C24[index]*dVy_dz +
                          C25[index]*dVz_dx + C24[index]*dVz_dy + C23[index]*dVz_dz);                    
    
        Tzz[index] += dt*(C13[index]*dVx_dx + C36[index]*dVx_dy + C35[index]*dVx_dz +
                          C36[index]*dVy_dx + C23[index]*dVy_dy + C34[index]*dVy_dz +
                          C35[index]*dVz_dx + C34[index]*dVz_dy + C33[index]*dVz_dz);  
    
        Txy[index] += dt*(C16[index]*dVx_dx + C66[index]*dVx_dy + C56[index]*dVx_dz +
                          C66[index]*dVy_dx + C26[index]*dVy_dy + C46[index]*dVy_dz +
                          C56[index]*dVz_dx + C46[index]*dVz_dy + C36[index]*dVz_dz);                    
    
        Txz[index] += dt*(C15[index]*dVx_dx + C56[index]*dVx_dy + C55[index]*dVx_dz +
                          C56[index]*dVy_dx + C25[index]*dVy_dy + C45[index]*dVy_dz +
                          C55[index]*dVz_dx + C45[index]*dVz_dy + C35[index]*dVz_dz);                    
    
        Tyz[index] += dt*(C14[index]*dVx_dx + C46[index]*dVx_dy + C45[index]*dVx_dz +
                          C46[index]*dVy_dx + C24[index]*dVy_dy + C44[index]*dVy_dz +
                          C45[index]*dVz_dx + C44[index]*dVz_dy + C34[index]*dVz_dz); 
    
        if ((i > 3) && (i < nzz-4) && (j > 3) && (j < nxx-4) && (k > 3) && (k < nyy-4))
        {
            P[index] = (Txx[index] + Tyy[index] + Tzz[index]) / 3.0f;
        }
    }
}

