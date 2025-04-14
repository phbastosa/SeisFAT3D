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

    auto * Cij = new float[nPoints]();

    std::string Cijkl_folder = catch_parameter("Cijkl_folder", parameters);

    auto * C11 = new float[volsize]();
    auto * uC11 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C11.bin", Cij, nPoints);
    expand_boundary(Cij, C11);
    compress(C11, uC11, nPoints, maxC11, minC11, compress_level);    
    cudaMalloc((void**)&(d_C11), volsize*sizeof(uintc));
    cudaMemcpy(d_C11, uC11, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C11;
    delete[] uC11;

    auto * C12 = new float[volsize]();
    auto * uC12 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C12.bin", Cij, nPoints);
    expand_boundary(Cij, C12);
    compress(C12, uC12, nPoints, maxC12, minC12, compress_level);    
    cudaMalloc((void**)&(d_C12), volsize*sizeof(uintc));
    cudaMemcpy(d_C12, uC12, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C12;
    delete[] uC12;

    auto * C13 = new float[volsize]();
    auto * uC13 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C13.bin", Cij, nPoints);
    expand_boundary(Cij, C13);
    compress(C13, uC13, nPoints, maxC13, minC13, compress_level);    
    cudaMalloc((void**)&(d_C13), volsize*sizeof(uintc));
    cudaMemcpy(d_C13, uC13, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C13;
    delete[] uC13;

    auto * C14 = new float[volsize]();
    auto * uC14 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C14.bin", Cij, nPoints);
    expand_boundary(Cij, C14);
    compress(C14, uC14, nPoints, maxC14, minC14, compress_level);    
    cudaMalloc((void**)&(d_C14), volsize*sizeof(uintc));
    cudaMemcpy(d_C14, uC14, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C14;
    delete[] uC14;

    auto * C15 = new float[volsize]();
    auto * uC15 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C15.bin", Cij, nPoints);
    expand_boundary(Cij, C15);
    compress(C15, uC15, nPoints, maxC15, minC15, compress_level);    
    cudaMalloc((void**)&(d_C15), volsize*sizeof(uintc));
    cudaMemcpy(d_C15, uC15, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C15;
    delete[] uC15;

    auto * C16 = new float[volsize]();
    auto * uC16 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C16.bin", Cij, nPoints);
    expand_boundary(Cij, C16);
    compress(C16, uC16, nPoints, maxC16, minC16, compress_level);    
    cudaMalloc((void**)&(d_C16), volsize*sizeof(uintc));
    cudaMemcpy(d_C16, uC16, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C16;
    delete[] uC16;

    auto * C22 = new float[volsize]();
    auto * uC22 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C22.bin", Cij, nPoints);
    expand_boundary(Cij, C22);
    compress(C22, uC22, nPoints, maxC22, minC22, compress_level);    
    cudaMalloc((void**)&(d_C22), volsize*sizeof(uintc));
    cudaMemcpy(d_C22, uC22, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C22;
    delete[] uC22;

    auto * C23 = new float[volsize]();
    auto * uC23 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C23.bin", Cij, nPoints);
    expand_boundary(Cij, C23);
    compress(C23, uC23, nPoints, maxC23, minC23, compress_level);    
    cudaMalloc((void**)&(d_C23), volsize*sizeof(uintc));
    cudaMemcpy(d_C23, uC23, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C23;
    delete[] uC23;
    
    auto * C24 = new float[volsize]();
    auto * uC24 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C24.bin", Cij, nPoints);
    expand_boundary(Cij, C24);
    compress(C24, uC24, nPoints, maxC24, minC24, compress_level);    
    cudaMalloc((void**)&(d_C24), volsize*sizeof(uintc));
    cudaMemcpy(d_C24, uC24, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C24;
    delete[] uC24;

    auto * C25 = new float[volsize]();
    auto * uC25 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C25.bin", Cij, nPoints);
    expand_boundary(Cij, C25);
    compress(C25, uC25, nPoints, maxC25, minC25, compress_level);    
    cudaMalloc((void**)&(d_C25), volsize*sizeof(uintc));
    cudaMemcpy(d_C25, uC25, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C25;
    delete[] uC25;

    auto * C26 = new float[volsize]();
    auto * uC26 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C26.bin", Cij, nPoints);
    expand_boundary(Cij, C26);
    compress(C26, uC26, nPoints, maxC26, minC26, compress_level);    
    cudaMalloc((void**)&(d_C26), volsize*sizeof(uintc));
    cudaMemcpy(d_C26, uC26, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C26;
    delete[] uC26;

    auto * C33 = new float[volsize]();
    auto * uC33 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C33.bin", Cij, nPoints);
    expand_boundary(Cij, C33);
    compress(C33, uC33, nPoints, maxC33, minC33, compress_level);    
    cudaMalloc((void**)&(d_C33), volsize*sizeof(uintc));
    cudaMemcpy(d_C33, uC33, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C33;
    delete[] uC33;
    
    auto * C34 = new float[volsize]();
    auto * uC34 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C34.bin", Cij, nPoints);
    expand_boundary(Cij, C34);
    compress(C34, uC34, nPoints, maxC34, minC34, compress_level);    
    cudaMalloc((void**)&(d_C34), volsize*sizeof(uintc));
    cudaMemcpy(d_C34, uC34, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C34;
    delete[] uC34;

    auto * C35 = new float[volsize]();
    auto * uC35 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C35.bin", Cij, nPoints);
    expand_boundary(Cij, C35);
    compress(C35, uC35, nPoints, maxC35, minC35, compress_level);    
    cudaMalloc((void**)&(d_C35), volsize*sizeof(uintc));
    cudaMemcpy(d_C35, uC35, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C35;
    delete[] uC35;

    auto * C36 = new float[volsize]();
    auto * uC36 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C36.bin", Cij, nPoints);
    expand_boundary(Cij, C36);
    compress(C36, uC36, nPoints, maxC36, minC36, compress_level);    
    cudaMalloc((void**)&(d_C36), volsize*sizeof(uintc));
    cudaMemcpy(d_C36, uC36, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C36;
    delete[] uC36;

    auto * C44 = new float[volsize]();
    auto * uC44 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C44.bin", Cij, nPoints);
    expand_boundary(Cij, C44);
    compress(C44, uC44, nPoints, maxC44, minC44, compress_level);    
    cudaMalloc((void**)&(d_C44), volsize*sizeof(uintc));
    cudaMemcpy(d_C44, uC44, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C44;
    delete[] uC44;

    auto * C45 = new float[volsize]();
    auto * uC45 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C45.bin", Cij, nPoints);
    expand_boundary(Cij, C45);
    compress(C45, uC45, nPoints, maxC45, minC45, compress_level);    
    cudaMalloc((void**)&(d_C45), volsize*sizeof(uintc));
    cudaMemcpy(d_C45, uC45, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C45;
    delete[] uC45;

    auto * C46 = new float[volsize]();
    auto * uC46 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C46.bin", Cij, nPoints);
    expand_boundary(Cij, C46);
    compress(C46, uC46, nPoints, maxC46, minC46, compress_level);    
    cudaMalloc((void**)&(d_C46), volsize*sizeof(uintc));
    cudaMemcpy(d_C46, uC46, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C46;
    delete[] uC46;

    auto * C55 = new float[volsize]();
    auto * uC55 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C55.bin", Cij, nPoints);
    expand_boundary(Cij, C55);
    compress(C55, uC55, nPoints, maxC55, minC55, compress_level);    
    cudaMalloc((void**)&(d_C55), volsize*sizeof(uintc));
    cudaMemcpy(d_C55, uC55, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C55;
    delete[] uC55;

    auto * C56 = new float[volsize]();
    auto * uC56 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C56.bin", Cij, nPoints);
    expand_boundary(Cij, C56);
    compress(C56, uC56, nPoints, maxC56, minC56, compress_level);    
    cudaMalloc((void**)&(d_C56), volsize*sizeof(uintc));
    cudaMemcpy(d_C56, uC56, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C56;
    delete[] uC56;

    auto * C66 = new float[volsize]();
    auto * uC66 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C66.bin", Cij, nPoints);
    expand_boundary(Cij, C66);
    compress(C66, uC66, nPoints, maxC66, minC66, compress_level);    
    cudaMalloc((void**)&(d_C66), volsize*sizeof(uintc));
    cudaMemcpy(d_C66, uC66, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C66;
    delete[] uC66;

    delete[] Cij;
}

void Elastic_ANI::propagation()
{
    for (int tId = 0; tId < nt + tlag; tId++)
    {
<<<<<<< HEAD
        compute_velocity_rsg<<<nBlocks, nThreads>>>(d_Vx, d_Vy, d_Vz, d_Txx, d_Tyy, d_Tzz, d_Txz, d_Tyz, d_Txy, d_B, d_T, 
                                                    d1D, d2D, d3D, wavelet, dwc, dx, dy, dz, dt, tId, tlag, sIdx, sIdy, 
                                                    sIdz, nxx, nyy, nzz, nb, nt);

        compute_pressure_rsg<<<nBlocks, nThreads>>>(d_Vx, d_Vy, d_Vz, d_Txx, d_Tyy, d_Tzz, d_Txz, d_Tyz, d_Txy, d_P, d_T, 
                                                    d_C11, d_C12, d_C13, d_C14, d_C15, d_C16, d_C22, d_C23, d_C24, d_C25, 
                                                    d_C26, d_C33, d_C34, d_C35, d_C36, d_C44, d_C45, d_C46, d_C55, d_C56, 
                                                    d_C66, tId, tlag, dx, dy, dz, dt, nxx, nyy, nzz, compress_level, 
                                                    minC11, maxC11, minC12, maxC12, minC13, maxC13, minC14, maxC14, 
                                                    minC15, maxC15, minC16, maxC16, minC22, maxC22, minC23, maxC23, 
                                                    minC24, maxC24, minC25, maxC25, minC26, maxC26, minC33, maxC33, 
                                                    minC34, maxC34, minC35, maxC35, minC36, maxC36, minC44, maxC44, 
                                                    minC45, maxC45, minC46, maxC46, minC55, maxC55, minC56, maxC56,
                                                    minC66, maxC66);
=======
        compute_velocity_rsg<<<nBlocks, nThreads>>>(d_Vx, d_Vy, d_Vz, d_Txx, d_Tyy, d_Tzz, d_Txz, d_Tyz, d_Txy, d_B, d_T, d1D, d2D, d3D, wavelet, dwc, dx, dy, dz, dt, tId, tlag, sIdx, sIdy, sIdz, nxx, nyy, nzz, nb, nt);

        compute_pressure_rsg<<<nBlocks, nThreads>>>(d_Vx, d_Vy, d_Vz, d_Txx, d_Tyy, d_Tzz, d_Txz, d_Tyz, d_Txy, d_P, d_T, d_C11, d_C12, d_C13, d_C14, d_C15, d_C16, d_C22, d_C23, d_C24, d_C25, d_C26, d_C33, d_C34, d_C35, d_C36, d_C44, d_C45, d_C46, d_C55, d_C56, d_C66, tId, tlag, dx, dy, dz, dt, nxx, nyy, nzz);
>>>>>>> b0e1613577e048fb2ad1900b573c24e644f53e66

        compute_seismogram<<<sBlocks, nThreads>>>(d_P, rIdx, rIdy, rIdz, seismogram, geometry->spread[srcId], tId, tlag, nt, nxx, nzz);     
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
                                     uintc * C11, uintc * C12, uintc * C13, uintc * C14, uintc * C15, uintc * C16, uintc * C22, uintc * C23, uintc * C24, uintc * C25, 
                                     uintc * C26, uintc * C33, uintc * C34, uintc * C35, uintc * C36, uintc * C44, uintc * C45, uintc * C46, uintc * C55, uintc * C56, 
                                     uintc * C66, int tId, int tlag, float dx, float dy, float dz, float dt, int nxx, int nyy, int nzz, int compress_level, float minC11, float maxC11, 
                                     float minC12, float maxC12, float minC13, float maxC13, float minC14, float maxC14, float minC15, float maxC15, float minC16, float maxC16, float minC22, 
                                     float maxC22, float minC23, float maxC23, float minC24, float maxC24, float minC25, float maxC25, float minC26, float maxC26, float minC33, float maxC33, 
                                     float minC34, float maxC34, float minC35, float maxC35, float minC36, float maxC36, float minC44, float maxC44, float minC45, float maxC45, float minC46, 
                                     float maxC46, float minC55, float maxC55, float minC56, float maxC56, float minC66, float maxC66)
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

        float c11 = minC11 + ((float)(C12[index])*(maxC11 - minC11) - 1) / (compress_level - 1);
        float c12 = minC12 + ((float)(C12[index])*(maxC12 - minC12) - 1) / (compress_level - 1);
        float c13 = minC13 + ((float)(C13[index])*(maxC13 - minC13) - 1) / (compress_level - 1);
        float c14 = minC14 + ((float)(C14[index])*(maxC14 - minC14) - 1) / (compress_level - 1);
        float c15 = minC15 + ((float)(C15[index])*(maxC15 - minC15) - 1) / (compress_level - 1);
        float c16 = minC16 + ((float)(C16[index])*(maxC16 - minC16) - 1) / (compress_level - 1);
    
        float c22 = minC22 + ((float)(C22[index])*(maxC22 - minC22) - 1) / (compress_level - 1);
        float c23 = minC23 + ((float)(C23[index])*(maxC23 - minC23) - 1) / (compress_level - 1);
        float c24 = minC24 + ((float)(C24[index])*(maxC24 - minC24) - 1) / (compress_level - 1);
        float c25 = minC25 + ((float)(C25[index])*(maxC25 - minC25) - 1) / (compress_level - 1);
        float c26 = minC26 + ((float)(C26[index])*(maxC26 - minC26) - 1) / (compress_level - 1);
    
        float c33 = minC33 + ((float)(C33[index])*(maxC33 - minC33) - 1) / (compress_level - 1);
        float c34 = minC34 + ((float)(C34[index])*(maxC34 - minC34) - 1) / (compress_level - 1);
        float c35 = minC35 + ((float)(C35[index])*(maxC35 - minC35) - 1) / (compress_level - 1);
        float c36 = minC36 + ((float)(C36[index])*(maxC36 - minC36) - 1) / (compress_level - 1);
    
        float c44 = minC44 + ((float)(C44[index])*(maxC44 - minC44) - 1) / (compress_level - 1);
        float c45 = minC45 + ((float)(C45[index])*(maxC45 - minC45) - 1) / (compress_level - 1);
        float c46 = minC46 + ((float)(C46[index])*(maxC46 - minC46) - 1) / (compress_level - 1);
    
        float c55 = minC55 + ((float)(C55[index])*(maxC55 - minC55) - 1) / (compress_level - 1);
        float c56 = minC56 + ((float)(C56[index])*(maxC56 - minC56) - 1) / (compress_level - 1);
    
        float c66 = minC66 + ((float)(C66[index])*(maxC66 - minC66) - 1) / (compress_level - 1);
            
        Txx[index] += dt*(c11*dVx_dx + c16*dVx_dy + c15*dVx_dz +
                          c16*dVy_dx + c12*dVy_dy + c14*dVy_dz +
                          c15*dVz_dx + c14*dVz_dy + c13*dVz_dz);                    
    
        Tyy[index] += dt*(c12*dVx_dx + c26*dVx_dy + c25*dVx_dz +
                          c26*dVy_dx + c22*dVy_dy + c24*dVy_dz +
                          c25*dVz_dx + c24*dVz_dy + c23*dVz_dz);                    
    
        Tzz[index] += dt*(c13*dVx_dx + c36*dVx_dy + c35*dVx_dz +
                          c36*dVy_dx + c23*dVy_dy + c34*dVy_dz +
                          c35*dVz_dx + c34*dVz_dy + c33*dVz_dz);  
    
        Txy[index] += dt*(c16*dVx_dx + c66*dVx_dy + c56*dVx_dz +
                          c66*dVy_dx + c26*dVy_dy + c46*dVy_dz +
                          c56*dVz_dx + c46*dVz_dy + c36*dVz_dz);                    
    
        Txz[index] += dt*(c15*dVx_dx + c56*dVx_dy + c55*dVx_dz +
                          c56*dVy_dx + c25*dVy_dy + c45*dVy_dz +
                          c55*dVz_dx + c45*dVz_dy + c35*dVz_dz);                    
    
        Tyz[index] += dt*(c14*dVx_dx + c46*dVx_dy + c45*dVx_dz +
                          c46*dVy_dx + c24*dVy_dy + c44*dVy_dz +
                          c45*dVz_dx + c44*dVz_dy + c34*dVz_dz); 
    
        if ((i > 3) && (i < nzz-4) && (j > 3) && (j < nxx-4) && (k > 3) && (k < nyy-4))
        {
            P[index] = (Txx[index] + Tyy[index] + Tzz[index]) / 3.0f;
        }
    }
}

