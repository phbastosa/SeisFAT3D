# include "eikonal_ani.cuh"

void Eikonal_ANI::set_properties()
{
    Vp = new float[nPoints]();
    Ro = new float[nPoints]();

    S = new float[volsize]();
    qS = new float[volsize]();

    std::string vp_file = catch_parameter("vp_model_file", parameters);
    std::string ro_file = catch_parameter("ro_model_file", parameters);

    import_binary_float(vp_file, Vp, nPoints);
    import_binary_float(ro_file, Ro, nPoints);

    float * slowness = new float[nPoints]();

    # pragma omp parallel for
    for (int index = 0; index < nPoints; index++)
        slowness[index] = 1.0f / Vp[index];

    expand_boundary(slowness, S);

    # pragma omp parallel for
    for (int index = 0; index < volsize; index++)
        qS[index] = S[index];

    delete[] slowness;
}

void Eikonal_ANI::set_conditions()
{
    modeling_type = "eikonal_ani";
    modeling_name = "Modeling type: Anisotropic eikonal solver";

    n = 3;
    v = 6;

    p = new float[n]();
    C = new float[v*v]();
    G = new float[n*n]();
    Gv = new float[n]();

    float * Cij = new float[nPoints]();

    std::string Cijkl_folder = catch_parameter("Cijkl_folder", parameters);

    C11 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C11.bin", Cij, nPoints);
    compress(Cij, C11, nPoints, maxC11, minC11, compress_level);

    C12 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C12.bin", Cij, nPoints);
    compress(Cij, C12, nPoints, maxC12, minC12, compress_level);

    C13 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C13.bin", Cij, nPoints);
    compress(Cij, C13, nPoints, maxC13, minC13, compress_level);

    C14 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C14.bin", Cij, nPoints);
    compress(Cij, C14, nPoints, maxC14, minC14, compress_level);

    C15 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C15.bin", Cij, nPoints);
    compress(Cij, C15, nPoints, maxC15, minC15, compress_level);

    C16 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C16.bin", Cij, nPoints);
    compress(Cij, C16, nPoints, maxC16, minC16, compress_level);

    C22 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C22.bin", Cij, nPoints);
    compress(Cij, C22, nPoints, maxC22, minC22, compress_level);

    C23 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C23.bin", Cij, nPoints);
    compress(Cij, C23, nPoints, maxC23, minC23, compress_level);
    
    C24 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C24.bin", Cij, nPoints);
    compress(Cij, C24, nPoints, maxC24, minC24, compress_level);

    C25 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C25.bin", Cij, nPoints);
    compress(Cij, C25, nPoints, maxC25, minC25, compress_level);

    C26 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C26.bin", Cij, nPoints);
    compress(Cij, C26, nPoints, maxC26, minC26, compress_level);

    C33 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C33.bin", Cij, nPoints);
    compress(Cij, C33, nPoints, maxC33, minC33, compress_level);
    
    C34 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C34.bin", Cij, nPoints);
    compress(Cij, C34, nPoints, maxC34, minC34, compress_level);

    C35 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C35.bin", Cij, nPoints);
    compress(Cij, C35, nPoints, maxC35, minC35, compress_level);

    C36 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C36.bin", Cij, nPoints);
    compress(Cij, C36, nPoints, maxC36, minC36, compress_level);

    C44 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C44.bin", Cij, nPoints);
    compress(Cij, C44, nPoints, maxC44, minC44, compress_level);

    C45 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C45.bin", Cij, nPoints);
    compress(Cij, C45, nPoints, maxC45, minC45, compress_level);

    C46 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C46.bin", Cij, nPoints);
    compress(Cij, C46, nPoints, maxC46, minC46, compress_level);

    C55 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C55.bin", Cij, nPoints);
    compress(Cij, C55, nPoints, maxC55, minC55, compress_level);

    C56 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C56.bin", Cij, nPoints);
    compress(Cij, C56, nPoints, maxC56, minC56, compress_level);

    C66 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C66.bin", Cij, nPoints);
    compress(Cij, C66, nPoints, maxC66, minC66, compress_level);
}

void Eikonal_ANI::forward_solver()
{
    initialization();
    
    propagation();

    for (int index = 0; index < volsize; index++)
    {
        int k = (int) (index / (nxx*nzz));         
        int j = (int) (index - k*nxx*nzz) / nzz;    
        int i = (int) (index - j*nzz - k*nxx*nzz);  

        if ((i == sIdz) && (j == sIdx) && (k == sIdy))    
            continue;

        if ((i >= nb) && (i < nzz-nb) && (j >= nb) && (j < nxx-nb) && (k >= nb) && (k < nyy-nb))
        {
            aId = (i - nb) + (j - nb)*nz + (k - nb)*nx*nz;

            float dTz = 0.5f*(T[(i+1) + j*nzz + k*nxx*nzz] - T[(i-1) + j*nzz + k*nxx*nzz]) / dz;
            float dTx = 0.5f*(T[i + (j+1)*nzz + k*nxx*nzz] - T[i + (j-1)*nzz + k*nxx*nzz]) / dx;
            float dTy = 0.5f*(T[i + j*nzz + (k+1)*nxx*nzz] - T[i + j*nzz + (k-1)*nxx*nzz]) / dy;

            float norm = sqrtf(dTx*dTx + dTy*dTy + dTz*dTz);

            p[0] = dTx / norm;
            p[1] = dTy / norm;
            p[2] = dTz / norm;

            get_stiffness();
            get_christoffel();
            get_eigen_values();

            S[index] = 1.0f / sqrtf(Gv[0] * Ro[aId]);
        }
    }

    initialization();

    propagation();

    compute_seismogram();

    # pragma omp parallel for
    for (int index = 0; index < volsize; index++)
        S[index] = qS[index];
}

void Eikonal_ANI::get_stiffness()
{
    float c11 = minC11 + ((float)(C12[aId])*(maxC11 - minC11) - 1) / (compress_level - 1);
    float c12 = minC12 + ((float)(C12[aId])*(maxC12 - minC12) - 1) / (compress_level - 1);
    float c13 = minC13 + ((float)(C13[aId])*(maxC13 - minC13) - 1) / (compress_level - 1);
    float c14 = minC14 + ((float)(C14[aId])*(maxC14 - minC14) - 1) / (compress_level - 1);
    float c15 = minC15 + ((float)(C15[aId])*(maxC15 - minC15) - 1) / (compress_level - 1);
    float c16 = minC16 + ((float)(C16[aId])*(maxC16 - minC16) - 1) / (compress_level - 1);

    float c22 = minC22 + ((float)(C22[aId])*(maxC22 - minC22) - 1) / (compress_level - 1);
    float c23 = minC23 + ((float)(C23[aId])*(maxC23 - minC23) - 1) / (compress_level - 1);
    float c24 = minC24 + ((float)(C24[aId])*(maxC24 - minC24) - 1) / (compress_level - 1);
    float c25 = minC25 + ((float)(C25[aId])*(maxC25 - minC25) - 1) / (compress_level - 1);
    float c26 = minC26 + ((float)(C26[aId])*(maxC26 - minC26) - 1) / (compress_level - 1);

    float c33 = minC33 + ((float)(C33[aId])*(maxC33 - minC33) - 1) / (compress_level - 1);
    float c34 = minC34 + ((float)(C34[aId])*(maxC34 - minC34) - 1) / (compress_level - 1);
    float c35 = minC35 + ((float)(C35[aId])*(maxC35 - minC35) - 1) / (compress_level - 1);
    float c36 = minC36 + ((float)(C36[aId])*(maxC36 - minC36) - 1) / (compress_level - 1);

    float c44 = minC44 + ((float)(C44[aId])*(maxC44 - minC44) - 1) / (compress_level - 1);
    float c45 = minC45 + ((float)(C45[aId])*(maxC45 - minC45) - 1) / (compress_level - 1);
    float c46 = minC46 + ((float)(C46[aId])*(maxC46 - minC46) - 1) / (compress_level - 1);

    float c55 = minC55 + ((float)(C55[aId])*(maxC55 - minC55) - 1) / (compress_level - 1);
    float c56 = minC56 + ((float)(C56[aId])*(maxC56 - minC56) - 1) / (compress_level - 1);

    float c66 = minC66 + ((float)(C66[aId])*(maxC66 - minC66) - 1) / (compress_level - 1);

    C[0+0*v] = c11; C[0+1*v] = c12; C[0+2*v] = c13; C[0+3*v] = c14; C[0+4*v] = c15; C[0+5*v] = c16;
    C[1+0*v] = c12; C[1+1*v] = c22; C[1+2*v] = c23; C[1+3*v] = c24; C[1+4*v] = c25; C[1+5*v] = c26;
    C[2+0*v] = c13; C[2+1*v] = c23; C[2+2*v] = c33; C[2+3*v] = c34; C[2+4*v] = c35; C[2+5*v] = c36;
    C[3+0*v] = c14; C[3+1*v] = c24; C[3+2*v] = c34; C[3+3*v] = c44; C[3+4*v] = c45; C[3+5*v] = c46;
    C[4+0*v] = c15; C[4+1*v] = c25; C[4+2*v] = c35; C[4+3*v] = c45; C[4+4*v] = c55; C[4+5*v] = c56;
    C[5+0*v] = c16; C[5+1*v] = c26; C[5+2*v] = c36; C[5+3*v] = c46; C[5+4*v] = c56; C[5+5*v] = c66;

    for (int i = 0; i < v*v; i++) C[i] *= 1.0f / Ro[aId] / Ro[aId];         
}

void Eikonal_ANI::get_christoffel()
{
    for (int index = 0; index < n*n; index++) 
        G[index] = 0.0f; 

    for (int i = 0; i < n; i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            for (int k = 0; k < n; k++) 
            {
                for (int l = 0; l < n; l++) 
                {
                    int I = voigt_map(i, k);
                    int J = voigt_map(j, l);

                    G[i + j*n] += C[I + J*v]*p[k]*p[l];
                }
            }
        }
    }
}

void Eikonal_ANI::get_eigen_values()
{
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
        float u = cbrt(-0.5f*q);

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
    
    if (Gv[0] < Gv[1]) std::swap(Gv[0],Gv[1]);
    if (Gv[1] < Gv[2]) std::swap(Gv[1],Gv[2]);
    if (Gv[0] < Gv[1]) std::swap(Gv[0],Gv[1]);    
}

int Eikonal_ANI::voigt_map(int i, int j)
{
    if (i == j)
        return i;
    if (((i == 1) && (j == 2)) || ((i == 2) && (j == 1)))
        return 3;
    if (((i == 2) && (j == 0)) || ((i == 0) && (j == 2)))
        return 4;
    if (((i == 0) && (j == 1)) || ((i == 1) && (j == 0)))
        return 5;
    
    return -1;
}
