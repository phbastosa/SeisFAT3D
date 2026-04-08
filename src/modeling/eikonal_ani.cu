# include "eikonal_ani.cuh"

void Eikonal_ANI::set_conditions()
{
    modeling_type = "eikonal_ani";
    modeling_name = "Modeling type: Anisotropic eikonal solver";

    auto * iModel = new float[nPoints]();
    auto * xModel = new float[volsize]();
    auto * uModel = new uintc[volsize]();

    import_binary_float(Cijkl_folder + "C11.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC11, minC11);        
    cudaMalloc((void**)&(d_C11), volsize*sizeof(uintc));
    cudaMemcpy(d_C11, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C12.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC12, minC12);        
    cudaMalloc((void**)&(d_C12), volsize*sizeof(uintc));
    cudaMemcpy(d_C12, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C13.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC13, minC13);        
    cudaMalloc((void**)&(d_C13), volsize*sizeof(uintc));
    cudaMemcpy(d_C13, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C14.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC14, minC14);        
    cudaMalloc((void**)&(d_C14), volsize*sizeof(uintc));
    cudaMemcpy(d_C14, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C15.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC15, minC15);        
    cudaMalloc((void**)&(d_C15), volsize*sizeof(uintc));
    cudaMemcpy(d_C15, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C16.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC16, minC16);        
    cudaMalloc((void**)&(d_C16), volsize*sizeof(uintc));
    cudaMemcpy(d_C16, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C22.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC22, minC22);        
    cudaMalloc((void**)&(d_C22), volsize*sizeof(uintc));
    cudaMemcpy(d_C22, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C23.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC23, minC23);        
    cudaMalloc((void**)&(d_C23), volsize*sizeof(uintc));
    cudaMemcpy(d_C23, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C24.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC24, minC24);        
    cudaMalloc((void**)&(d_C24), volsize*sizeof(uintc));
    cudaMemcpy(d_C24, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C25.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC25, minC25);        
    cudaMalloc((void**)&(d_C25), volsize*sizeof(uintc));
    cudaMemcpy(d_C25, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C26.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC26, minC26);        
    cudaMalloc((void**)&(d_C26), volsize*sizeof(uintc));
    cudaMemcpy(d_C26, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C33.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC33, minC33);        
    cudaMalloc((void**)&(d_C33), volsize*sizeof(uintc));
    cudaMemcpy(d_C33, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C34.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC34, minC34);        
    cudaMalloc((void**)&(d_C34), volsize*sizeof(uintc));
    cudaMemcpy(d_C34, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C35.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC35, minC35);        
    cudaMalloc((void**)&(d_C35), volsize*sizeof(uintc));
    cudaMemcpy(d_C35, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C36.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC36, minC36);        
    cudaMalloc((void**)&(d_C36), volsize*sizeof(uintc));
    cudaMemcpy(d_C36, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C44.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC44, minC44);        
    cudaMalloc((void**)&(d_C44), volsize*sizeof(uintc));
    cudaMemcpy(d_C44, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C45.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC45, minC45);        
    cudaMalloc((void**)&(d_C45), volsize*sizeof(uintc));
    cudaMemcpy(d_C45, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C46.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC46, minC46);        
    cudaMalloc((void**)&(d_C46), volsize*sizeof(uintc));
    cudaMemcpy(d_C46, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C55.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC55, minC55);        
    cudaMalloc((void**)&(d_C55), volsize*sizeof(uintc));
    cudaMemcpy(d_C55, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C56.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC56, minC56);        
    cudaMalloc((void**)&(d_C56), volsize*sizeof(uintc));
    cudaMemcpy(d_C56, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);

    import_binary_float(Cijkl_folder + "C66.bin", iModel, nPoints);
    expand_boundary(iModel, xModel);
    compression(xModel, uModel, volsize, maxC66, minC66);        
    cudaMalloc((void**)&(d_C66), volsize*sizeof(uintc));
    cudaMemcpy(d_C66, uModel, volsize*sizeof(uintc), cudaMemcpyHostToDevice);
    
    delete[] iModel;   
    delete[] xModel;
    delete[] uModel;
}

void Eikonal_ANI::time_propagation()
{
    initialization();
    eikonal_solver();

    get_quasi_slowness<<<nBlocks,NTHREADS>>>(d_T,d_S,dx,dy,dz,sIdx,sIdy,sIdz,nxx,nyy,nzz,nb,d_C11, 
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
    const float EPS = 1e-12f;

    const int n = 3;
    const int v = 6;

    float p[n];
    float C[v*v];
    float G[n*n];

    int index = blockIdx.x * blockDim.x + threadIdx.x;

    int k = index / (nxx * nzz);
    int j = (index - k * nxx * nzz) / nzz;
    int i = index - j * nzz - k * nxx * nzz;

    if (i < nb || i >= nzz - nb ||
        j < nb || j >= nxx - nb ||
        k < nb || k >= nyy - nb)
        return;

    if (i == sIdz && j == sIdx && k == sIdy)
        return;

    float dTz = 0.5f*(T[(i+1) + j*nzz + k*nxx*nzz] - T[(i-1) + j*nzz + k*nxx*nzz]) / dz;
    float dTx = 0.5f*(T[i + (j+1)*nzz + k*nxx*nzz] - T[i + (j-1)*nzz + k*nxx*nzz]) / dx;
    float dTy = 0.5f*(T[i + j*nzz + (k+1)*nxx*nzz] - T[i + j*nzz + (k-1)*nxx*nzz]) / dy;

    float norm = sqrtf(dTx*dTx + dTy*dTy + dTz*dTz) + EPS;

    p[0] = dTx / norm;
    p[1] = dTy / norm;
    p[2] = dTz / norm;

    auto decode = [&](uintc v, float vmin, float vmax) {
        return vmin + (float(v) - 1.0f) * (vmax - vmin) / (COMPRESS - 1);
    };

    float c11 = decode(C11[index], minC11, maxC11); float c23 = decode(C23[index], minC23, maxC23); float c36 = decode(C36[index], minC36, maxC36); 
    float c12 = decode(C12[index], minC12, maxC12); float c24 = decode(C24[index], minC24, maxC24); float c44 = decode(C44[index], minC44, maxC44); 
    float c13 = decode(C13[index], minC13, maxC13); float c25 = decode(C25[index], minC25, maxC25); float c45 = decode(C45[index], minC45, maxC45); 
    float c14 = decode(C14[index], minC14, maxC14); float c26 = decode(C26[index], minC26, maxC26); float c46 = decode(C46[index], minC46, maxC46); 
    float c15 = decode(C15[index], minC15, maxC15); float c33 = decode(C33[index], minC33, maxC33); float c55 = decode(C55[index], minC55, maxC55); 
    float c16 = decode(C16[index], minC16, maxC16); float c34 = decode(C34[index], minC34, maxC34); float c56 = decode(C56[index], minC56, maxC56); 
    float c22 = decode(C22[index], minC22, maxC22); float c35 = decode(C35[index], minC35, maxC35); float c66 = decode(C66[index], minC66, maxC66); 
    
    float s_val = S[index];
    float Ro = c33*s_val*s_val;

    const int voigt_map[n][n] = {{0, 5, 4}, {5, 1, 3}, {4, 3, 2}};

    C[0+0*v] = c11; C[0+1*v] = c12; C[0+2*v] = c13; C[0+3*v] = c14; C[0+4*v] = c15; C[0+5*v] = c16;
    C[1+0*v] = c12; C[1+1*v] = c22; C[1+2*v] = c23; C[1+3*v] = c24; C[1+4*v] = c25; C[1+5*v] = c26;
    C[2+0*v] = c13; C[2+1*v] = c23; C[2+2*v] = c33; C[2+3*v] = c34; C[2+4*v] = c35; C[2+5*v] = c36;
    C[3+0*v] = c14; C[3+1*v] = c24; C[3+2*v] = c34; C[3+3*v] = c44; C[3+4*v] = c45; C[3+5*v] = c46;
    C[4+0*v] = c15; C[4+1*v] = c25; C[4+2*v] = c35; C[4+3*v] = c45; C[4+4*v] = c55; C[4+5*v] = c56;
    C[5+0*v] = c16; C[5+1*v] = c26; C[5+2*v] = c36; C[5+3*v] = c46; C[5+4*v] = c56; C[5+5*v] = c66;

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

    G[1] = G[3] = 0.5*(G[1] + G[3]);
    G[2] = G[6] = 0.5*(G[2] + G[6]);
    G[5] = G[7] = 0.5*(G[5] + G[7]);

    double a = -(G[0] + G[4] + G[8]);
    double b = G[0]*G[4] + G[4]*G[8] + G[0]*G[8]
             - G[1]*G[3] - G[2]*G[6] - G[5]*G[7];

    double c = -(G[0]*(G[4]*G[8] - G[5]*G[7])
               - G[1]*(G[3]*G[8] - G[5]*G[6])
               + G[2]*(G[3]*G[7] - G[4]*G[6]));

    double pc = b - a*a/3.0;
    double qc = 2.0*a*a*a/27.0 - a*b/3.0 + c;

    double det = 0.25*qc*qc + (pc*pc*pc)/27.0;

    double Gv[n] = {0.0, 0.0, 0.0};

    if (det > EPS)
    {
        double sqrt_det = sqrt(det);
        double u = cbrt(-0.5*qc + sqrt_det);
        double v = cbrt(-0.5*qc - sqrt_det);
        Gv[0] = u + v - a/3.0;
    }
    else
    {
        double r = sqrt(max(-pc*pc*pc/27.0, 0.0));
        double arg = -0.5*qc / max(r, EPS);
        arg = min(1.0, max(-1.0, arg));

        double phi = acos(arg);
        r = 2.0 * cbrt(r);

        Gv[0] = r*cos(phi/3.0) - a/3.0;
        Gv[1] = r*cos((phi+2.0*M_PI)/3.0) - a/3.0;
        Gv[2] = r*cos((phi+4.0*M_PI)/3.0) - a/3.0;
    }

    double lambda_max = max(Gv[0], max(Gv[1], Gv[2]));
    lambda_max = max(lambda_max, EPS);

    S[index] = 1.0f / sqrtf((float)(lambda_max * Ro));
}
