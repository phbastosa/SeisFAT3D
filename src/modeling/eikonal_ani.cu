# include "eikonal_ani.cuh"

void Eikonal_ANI::set_conditions()
{
    modeling_type = "eikonal_ani";
    modeling_name = "Modeling type: Anisotropic eikonal solver";

    auto * Cij = new float[nPoints]();

    std::string Cijkl_folder = catch_parameter("Cijkl_folder", parameters);

    auto * C11 = new float[volsize]();
    auto * uC11 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C11.bin", Cij, nPoints);
    expand_boundary(Cij, C11);
    compression(C11, uC11, volsize, maxC11, minC11);        
    cudaMalloc((void**)&(d_C11), volsize*sizeof(uintc));
    cudaMemcpy(d_C11, uC11, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C11;
    delete[] uC11;

    auto * C12 = new float[volsize]();
    auto * uC12 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C12.bin", Cij, nPoints);
    expand_boundary(Cij, C12);
    compression(C12, uC12, volsize, maxC12, minC12);        
    cudaMalloc((void**)&(d_C12), volsize*sizeof(uintc));
    cudaMemcpy(d_C12, uC12, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C12;
    delete[] uC12;

    auto * C13 = new float[volsize]();
    auto * uC13 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C13.bin", Cij, nPoints);
    expand_boundary(Cij, C13);
    compression(C13, uC13, volsize, maxC13, minC13);    
    cudaMalloc((void**)&(d_C13), volsize*sizeof(uintc));
    cudaMemcpy(d_C13, uC13, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C13;
    delete[] uC13;

    auto * C14 = new float[volsize]();
    auto * uC14 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C14.bin", Cij, nPoints);
    expand_boundary(Cij, C14);
    compression(C14, uC14, volsize, maxC14, minC14);    
    cudaMalloc((void**)&(d_C14), volsize*sizeof(uintc));
    cudaMemcpy(d_C14, uC14, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C14;
    delete[] uC14;

    auto * C15 = new float[volsize]();
    auto * uC15 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C15.bin", Cij, nPoints);
    expand_boundary(Cij, C15);
    compression(C15, uC15, volsize, maxC15, minC15);    
    cudaMalloc((void**)&(d_C15), volsize*sizeof(uintc));
    cudaMemcpy(d_C15, uC15, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C15;
    delete[] uC15;

    auto * C16 = new float[volsize]();
    auto * uC16 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C16.bin", Cij, nPoints);
    expand_boundary(Cij, C16);
    compression(C16, uC16, volsize, maxC16, minC16);    
    cudaMalloc((void**)&(d_C16), volsize*sizeof(uintc));
    cudaMemcpy(d_C16, uC16, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C16;
    delete[] uC16;

    auto * C22 = new float[volsize]();
    auto * uC22 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C22.bin", Cij, nPoints);
    expand_boundary(Cij, C22);
    compression(C22, uC22, volsize, maxC22, minC22);    
    cudaMalloc((void**)&(d_C22), volsize*sizeof(uintc));
    cudaMemcpy(d_C22, uC22, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C22;
    delete[] uC22;

    auto * C23 = new float[volsize]();
    auto * uC23 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C23.bin", Cij, nPoints);
    expand_boundary(Cij, C23);
    compression(C23, uC23, volsize, maxC23, minC23);    
    cudaMalloc((void**)&(d_C23), volsize*sizeof(uintc));
    cudaMemcpy(d_C23, uC23, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C23;
    delete[] uC23;
    
    auto * C24 = new float[volsize]();
    auto * uC24 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C24.bin", Cij, nPoints);
    expand_boundary(Cij, C24);
    compression(C24, uC24, volsize, maxC24, minC24);    
    cudaMalloc((void**)&(d_C24), volsize*sizeof(uintc));
    cudaMemcpy(d_C24, uC24, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C24;
    delete[] uC24;

    auto * C25 = new float[volsize]();
    auto * uC25 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C25.bin", Cij, nPoints);
    expand_boundary(Cij, C25);
    compression(C25, uC25, volsize, maxC25, minC25);    
    cudaMalloc((void**)&(d_C25), volsize*sizeof(uintc));
    cudaMemcpy(d_C25, uC25, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C25;
    delete[] uC25;

    auto * C26 = new float[volsize]();
    auto * uC26 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C26.bin", Cij, nPoints);
    expand_boundary(Cij, C26);
    compression(C26, uC26, volsize, maxC26, minC26);    
    cudaMalloc((void**)&(d_C26), volsize*sizeof(uintc));
    cudaMemcpy(d_C26, uC26, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C26;
    delete[] uC26;

    auto * C33 = new float[volsize]();
    auto * uC33 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C33.bin", Cij, nPoints);
    expand_boundary(Cij, C33);
    compression(C33, uC33, volsize, maxC33, minC33);    
    cudaMalloc((void**)&(d_C33), volsize*sizeof(uintc));
    cudaMemcpy(d_C33, uC33, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C33;
    delete[] uC33;
    
    auto * C34 = new float[volsize]();
    auto * uC34 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C34.bin", Cij, nPoints);
    expand_boundary(Cij, C34);
    compression(C34, uC34, volsize, maxC34, minC34);    
    cudaMalloc((void**)&(d_C34), volsize*sizeof(uintc));
    cudaMemcpy(d_C34, uC34, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C34;
    delete[] uC34;

    auto * C35 = new float[volsize]();
    auto * uC35 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C35.bin", Cij, nPoints);
    expand_boundary(Cij, C35);
    compression(C35, uC35, volsize, maxC35, minC35);    
    cudaMalloc((void**)&(d_C35), volsize*sizeof(uintc));
    cudaMemcpy(d_C35, uC35, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C35;
    delete[] uC35;

    auto * C36 = new float[volsize]();
    auto * uC36 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C36.bin", Cij, nPoints);
    expand_boundary(Cij, C36);
    compression(C36, uC36, volsize, maxC36, minC36);    
    cudaMalloc((void**)&(d_C36), volsize*sizeof(uintc));
    cudaMemcpy(d_C36, uC36, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C36;
    delete[] uC36;

    auto * C44 = new float[volsize]();
    auto * uC44 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C44.bin", Cij, nPoints);
    expand_boundary(Cij, C44);
    compression(C44, uC44, volsize, maxC44, minC44);    
    cudaMalloc((void**)&(d_C44), volsize*sizeof(uintc));
    cudaMemcpy(d_C44, uC44, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C44;
    delete[] uC44;

    auto * C45 = new float[volsize]();
    auto * uC45 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C45.bin", Cij, nPoints);
    expand_boundary(Cij, C45);
    compression(C45, uC45, volsize, maxC45, minC45);    
    cudaMalloc((void**)&(d_C45), volsize*sizeof(uintc));
    cudaMemcpy(d_C45, uC45, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C45;
    delete[] uC45;

    auto * C46 = new float[volsize]();
    auto * uC46 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C46.bin", Cij, nPoints);
    expand_boundary(Cij, C46);
    compression(C46, uC46, volsize, maxC46, minC46);    
    cudaMalloc((void**)&(d_C46), volsize*sizeof(uintc));
    cudaMemcpy(d_C46, uC46, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C46;
    delete[] uC46;

    auto * C55 = new float[volsize]();
    auto * uC55 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C55.bin", Cij, nPoints);
    expand_boundary(Cij, C55);
    compression(C55, uC55, volsize, maxC55, minC55);    
    cudaMalloc((void**)&(d_C55), volsize*sizeof(uintc));
    cudaMemcpy(d_C55, uC55, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C55;
    delete[] uC55;

    auto * C56 = new float[volsize]();
    auto * uC56 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C56.bin", Cij, nPoints);
    expand_boundary(Cij, C56);
    compression(C56, uC56, volsize, maxC56, minC56);    
    cudaMalloc((void**)&(d_C56), volsize*sizeof(uintc));
    cudaMemcpy(d_C56, uC56, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C56;
    delete[] uC56;

    auto * C66 = new float[volsize]();
    auto * uC66 = new uintc[volsize]();
    import_binary_float(Cijkl_folder + "C66.bin", Cij, nPoints);
    expand_boundary(Cij, C66);
    compression(C66, uC66, volsize, maxC66, minC66);    
    cudaMalloc((void**)&(d_C66), volsize*sizeof(uintc));
    cudaMemcpy(d_C66, uC66, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] C66;
    delete[] uC66;

    delete[] Cij;
}

void Eikonal_ANI::time_propagation()
{
    initialization();
    eikonal_solver();

    get_quasi_slowness<<<nBlocks,nThreads>>>(d_T,d_S,dx,dy,dz,sIdx,sIdy,sIdz,nxx,nyy,nzz,nb,d_C11, 
                                             d_C12,d_C13,d_C14,d_C15,d_C16,d_C22,d_C23,d_C24,d_C25, 
                                             d_C26,d_C33,d_C34,d_C35,d_C36,d_C44,d_C45,d_C46,d_C55, 
                                             d_C56,d_C66,minC11,maxC11,minC12,maxC12,minC13,maxC13,
                                             minC14,maxC14,minC15,maxC15,minC16,maxC16,minC22,maxC22,
                                             minC23,maxC23,minC24,maxC24,minC25,maxC25,minC26,maxC26,
                                             minC33,maxC33,minC34,maxC34,minC35,maxC35,minC36,maxC36,
                                             minC44,maxC44,minC45,maxC45,minC46,maxC46,minC55,maxC55,
                                             minC56,maxC56,minC66,maxC66);
                                             
    initialization();
    eikonal_solver();
    compute_seismogram();

    cudaMemcpy(d_S, S, volsize * sizeof(float), cudaMemcpyHostToDevice);
}

__global__ void get_quasi_slowness(float * T, float * S, float dx, float dy, float dz, int sIdx, int sIdy, int sIdz, int nxx, int nyy, int nzz, int nb,
                                   uintc * C11, uintc * C12, uintc * C13, uintc * C14, uintc * C15, uintc * C16, uintc * C22, uintc * C23, uintc * C24, uintc * C25, 
                                   uintc * C26, uintc * C33, uintc * C34, uintc * C35, uintc * C36, uintc * C44, uintc * C45, uintc * C46, uintc * C55, uintc * C56, 
                                   uintc * C66, float minC11, float maxC11, float minC12, float maxC12, float minC13, float maxC13, float minC14, float maxC14, 
                                   float minC15, float maxC15, float minC16, float maxC16, float minC22, float maxC22, float minC23, float maxC23, float minC24, 
                                   float maxC24, float minC25, float maxC25, float minC26, float maxC26, float minC33, float maxC33, float minC34, float maxC34, 
                                   float minC35, float maxC35, float minC36, float maxC36, float minC44, float maxC44, float minC45, float maxC45, float minC46, 
                                   float maxC46, float minC55, float maxC55, float minC56, float maxC56, float minC66, float maxC66)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = (int) (index / (nxx*nzz));         
    int j = (int) (index - k*nxx*nzz) / nzz;    
    int i = (int) (index - j*nzz - k*nxx*nzz);  

    const int n = 3;
    const int v = 6;

    float p[n];
    float C[v*v];
    float G[n*n];
    float Gv[n];

    int voigt_map[n][n] = {{0, 5, 4}, {5, 1, 3}, {4, 3, 2}};

    if ((i >= nb) && (i < nzz-nb) && (j >= nb) && (j < nxx-nb) && (k >= nb) && (k < nyy-nb))
    {
        if (!((i == sIdz) && (j == sIdx) && (k == sIdy)))    
        {
            float dTz = 0.5f*(T[(i+1) + j*nzz + k*nxx*nzz] - T[(i-1) + j*nzz + k*nxx*nzz]) / dz;
            float dTx = 0.5f*(T[i + (j+1)*nzz + k*nxx*nzz] - T[i + (j-1)*nzz + k*nxx*nzz]) / dx;
            float dTy = 0.5f*(T[i + j*nzz + (k+1)*nxx*nzz] - T[i + j*nzz + (k-1)*nxx*nzz]) / dy;

            float norm = sqrtf(dTx*dTx + dTy*dTy + dTz*dTz);

            p[0] = dTx / norm;
            p[1] = dTy / norm;
            p[2] = dTz / norm;
            
            float c11 = (minC11 + (static_cast<float>(C11[index]) - 1.0f) * (maxC11 - minC11) / (COMPRESS - 1));
            float c12 = (minC12 + (static_cast<float>(C12[index]) - 1.0f) * (maxC12 - minC12) / (COMPRESS - 1));
            float c13 = (minC13 + (static_cast<float>(C13[index]) - 1.0f) * (maxC13 - minC13) / (COMPRESS - 1));
            float c14 = (minC14 + (static_cast<float>(C14[index]) - 1.0f) * (maxC14 - minC14) / (COMPRESS - 1));
            float c15 = (minC15 + (static_cast<float>(C15[index]) - 1.0f) * (maxC15 - minC15) / (COMPRESS - 1));
            float c16 = (minC16 + (static_cast<float>(C16[index]) - 1.0f) * (maxC16 - minC16) / (COMPRESS - 1));

            float c22 = (minC22 + (static_cast<float>(C22[index]) - 1.0f) * (maxC22 - minC22) / (COMPRESS - 1));
            float c23 = (minC23 + (static_cast<float>(C23[index]) - 1.0f) * (maxC23 - minC23) / (COMPRESS - 1));
            float c24 = (minC24 + (static_cast<float>(C24[index]) - 1.0f) * (maxC24 - minC24) / (COMPRESS - 1));
            float c25 = (minC25 + (static_cast<float>(C25[index]) - 1.0f) * (maxC25 - minC25) / (COMPRESS - 1));
            float c26 = (minC26 + (static_cast<float>(C26[index]) - 1.0f) * (maxC26 - minC26) / (COMPRESS - 1));

            float c33 = (minC33 + (static_cast<float>(C33[index]) - 1.0f) * (maxC33 - minC33) / (COMPRESS - 1));
            float c34 = (minC34 + (static_cast<float>(C34[index]) - 1.0f) * (maxC34 - minC34) / (COMPRESS - 1));
            float c35 = (minC35 + (static_cast<float>(C35[index]) - 1.0f) * (maxC35 - minC35) / (COMPRESS - 1));
            float c36 = (minC36 + (static_cast<float>(C36[index]) - 1.0f) * (maxC36 - minC36) / (COMPRESS - 1));

            float c44 = (minC44 + (static_cast<float>(C44[index]) - 1.0f) * (maxC44 - minC44) / (COMPRESS - 1));
            float c45 = (minC45 + (static_cast<float>(C45[index]) - 1.0f) * (maxC45 - minC45) / (COMPRESS - 1));
            float c46 = (minC46 + (static_cast<float>(C46[index]) - 1.0f) * (maxC46 - minC46) / (COMPRESS - 1));

            float c55 = (minC55 + (static_cast<float>(C55[index]) - 1.0f) * (maxC55 - minC55) / (COMPRESS - 1));
            float c56 = (minC56 + (static_cast<float>(C56[index]) - 1.0f) * (maxC56 - minC56) / (COMPRESS - 1));

            float c66 = (minC66 + (static_cast<float>(C66[index]) - 1.0f) * (maxC66 - minC66) / (COMPRESS - 1));

            C[0+0*v] = c11; C[0+1*v] = c12; C[0+2*v] = c13; C[0+3*v] = c14; C[0+4*v] = c15; C[0+5*v] = c16;
            C[1+0*v] = c12; C[1+1*v] = c22; C[1+2*v] = c23; C[1+3*v] = c24; C[1+4*v] = c25; C[1+5*v] = c26;
            C[2+0*v] = c13; C[2+1*v] = c23; C[2+2*v] = c33; C[2+3*v] = c34; C[2+4*v] = c35; C[2+5*v] = c36;
            C[3+0*v] = c14; C[3+1*v] = c24; C[3+2*v] = c34; C[3+3*v] = c44; C[3+4*v] = c45; C[3+5*v] = c46;
            C[4+0*v] = c15; C[4+1*v] = c25; C[4+2*v] = c35; C[4+3*v] = c45; C[4+4*v] = c55; C[4+5*v] = c56;
            C[5+0*v] = c16; C[5+1*v] = c26; C[5+2*v] = c36; C[5+3*v] = c46; C[5+4*v] = c56; C[5+5*v] = c66;

            float Ro = c33*S[index]*S[index];    
            
            for (int indp = 0; indp < v*v; indp++)
                C[indp] = C[indp] / Ro / Ro;

            for (int indp = 0; indp < n*n; indp++) 
                G[indp] = 0.0f; 

            for (int ip = 0; ip < n; ip++) 
            {
                for (int jp = 0; jp < n; jp++) 
                {
                    for (int kp = 0; kp < n; kp++) 
                    {
                        for (int lp = 0; lp < n; lp++) 
                        {
                            int I = voigt_map[ip][kp];
                            int J = voigt_map[jp][lp];

                            G[ip + jp*n] += C[I + J*v]*p[kp]*p[lp];
                        }
                    }
                }
            }

            float a = -(G[0] + G[4] + G[8]);
    
            float b = G[0]*G[4] + G[4]*G[8] + 
                      G[0]*G[8] - G[3]*G[1] - 
                      G[6]*G[6] - G[7]*G[5];
            
            float c = -(G[0]*(G[4]*G[8] - G[7]*G[5]) -
                        G[3]*(G[1]*G[8] - G[7]*G[6]) +
                        G[6]*(G[1]*G[5] - G[4]*G[6]));

            float p = b - (a*a)/3.0f;
            float q = (2.0f*a*a*a)/27.0f - (a*b)/3.0f + c;

            float detG = 0.25f*(q*q) + (p*p*p)/27.0f;

            if (detG > 0) 
            {
                float u = cbrtf(-0.5f*q + sqrtf(detG));
                float v = cbrtf(-0.5f*q - sqrtf(detG));
                
                Gv[0] = u + v - a/3.0f;
            } 
            else if (detG == 0) 
            {       
                float u = cbrtf(-0.5f*q);

                Gv[0] = 2.0f*u - a/3.0f;
                Gv[1] =-1.0f*u - a/3.0f;         
            } 
            else  
            {
                float r = sqrtf(-p*p*p/27.0f);
                float phi = acosf(-0.5f*q/r);
                
                r = 2.0f*cbrtf(r);

                Gv[0] = r*cosf(phi/3.0f) - a/3.0f;
                Gv[1] = r*cosf((phi + 2.0f*M_PI)/3.0f) - a/3.0f;  
                Gv[2] = r*cosf((phi + 4.0f*M_PI)/3.0f) - a/3.0f;      
            }
            
            float aux;

            if (Gv[0] < Gv[1]) {aux = Gv[0]; Gv[0] = Gv[1]; Gv[1] = aux;} 
            if (Gv[1] < Gv[2]) {aux = Gv[1]; Gv[1] = Gv[2]; Gv[2] = aux;}
            if (Gv[0] < Gv[1]) {aux = Gv[0]; Gv[0] = Gv[1]; Gv[1] = aux;}    

            S[index] = 1.0f / sqrtf(Gv[0] * Ro);
        }
    }
}
