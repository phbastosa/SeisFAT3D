# include "eikonal_ani.cuh"

void Eikonal_ANI::set_stiffness_element(std::string element, uintc * dCij, float &max, float &min)
{   
    std::string Cijkl_folder = catch_parameter("Cijkl_folder", parameters);

    auto * Cij = new float[volsize]();
    auto * uCij = new uintc[volsize]();
    auto * Caux = new float[nPoints]();
    import_binary_float(Cijkl_folder + element + ".bin", Caux, nPoints);
    expand_boundary(Caux, Cij);
    compression(Cij, uCij, volsize, max, min);        
    cudaMalloc((void**)&(dCij), volsize*sizeof(uintc));
    cudaMemcpy(dCij, uCij, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    delete[] Cij;
    delete[] uCij;
    delete[] Caux;
}

void Eikonal_ANI::set_conditions()
{
    modeling_type = "eikonal_ani";
    modeling_name = "Modeling type: Anisotropic eikonal solver";

    set_stiffness_element("C11", d_C11, maxC11, minC11);
    set_stiffness_element("C12", d_C12, maxC12, minC12);
    set_stiffness_element("C13", d_C13, maxC13, minC13);
    set_stiffness_element("C14", d_C14, maxC14, minC14);
    set_stiffness_element("C15", d_C15, maxC15, minC15);
    set_stiffness_element("C16", d_C16, maxC16, minC16);

    set_stiffness_element("C22", d_C22, maxC22, minC22);
    set_stiffness_element("C23", d_C23, maxC23, minC23);
    set_stiffness_element("C24", d_C24, maxC24, minC24);
    set_stiffness_element("C25", d_C25, maxC25, minC25);
    set_stiffness_element("C26", d_C26, maxC26, minC26);

    set_stiffness_element("C33", d_C33, maxC33, minC33);
    set_stiffness_element("C34", d_C34, maxC34, minC34);
    set_stiffness_element("C35", d_C35, maxC35, minC35);
    set_stiffness_element("C36", d_C36, maxC36, minC36);

    set_stiffness_element("C44", d_C44, maxC44, minC44);
    set_stiffness_element("C45", d_C45, maxC45, minC45);
    set_stiffness_element("C46", d_C46, maxC46, minC46);

    set_stiffness_element("C55", d_C55, maxC55, minC55);
    set_stiffness_element("C56", d_C56, maxC56, minC56);

    set_stiffness_element("C66", d_C66, maxC66, minC66);
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

    copy_slowness_to_device();
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
