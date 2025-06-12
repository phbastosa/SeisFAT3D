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
    compression(Cij, C11, nPoints, maxC11, minC11);

    C12 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C12.bin", Cij, nPoints);
    compression(Cij, C12, nPoints, maxC12, minC12);

    C13 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C13.bin", Cij, nPoints);
    compression(Cij, C13, nPoints, maxC13, minC13);

    C14 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C14.bin", Cij, nPoints);
    compression(Cij, C14, nPoints, maxC14, minC14);

    C15 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C15.bin", Cij, nPoints);
    compression(Cij, C15, nPoints, maxC15, minC15);

    C16 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C16.bin", Cij, nPoints);
    compression(Cij, C16, nPoints, maxC16, minC16);

    C22 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C22.bin", Cij, nPoints);
    compression(Cij, C22, nPoints, maxC22, minC22);

    C23 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C23.bin", Cij, nPoints);
    compression(Cij, C23, nPoints, maxC23, minC23);
    
    C24 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C24.bin", Cij, nPoints);
    compression(Cij, C24, nPoints, maxC24, minC24);

    C25 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C25.bin", Cij, nPoints);
    compression(Cij, C25, nPoints, maxC25, minC25);

    C26 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C26.bin", Cij, nPoints);
    compression(Cij, C26, nPoints, maxC26, minC26);

    C33 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C33.bin", Cij, nPoints);
    compression(Cij, C33, nPoints, maxC33, minC33);
    
    C34 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C34.bin", Cij, nPoints);
    compression(Cij, C34, nPoints, maxC34, minC34);

    C35 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C35.bin", Cij, nPoints);
    compression(Cij, C35, nPoints, maxC35, minC35);

    C36 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C36.bin", Cij, nPoints);
    compression(Cij, C36, nPoints, maxC36, minC36);

    C44 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C44.bin", Cij, nPoints);
    compression(Cij, C44, nPoints, maxC44, minC44);

    C45 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C45.bin", Cij, nPoints);
    compression(Cij, C45, nPoints, maxC45, minC45);

    C46 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C46.bin", Cij, nPoints);
    compression(Cij, C46, nPoints, maxC46, minC46);

    C55 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C55.bin", Cij, nPoints);
    compression(Cij, C55, nPoints, maxC55, minC55);

    C56 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C56.bin", Cij, nPoints);
    compression(Cij, C56, nPoints, maxC56, minC56);

    C66 = new uintc[nPoints]();
    import_binary_float(Cijkl_folder + "C66.bin", Cij, nPoints);
    compression(Cij, C66, nPoints, maxC66, minC66);
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
    float c11 = (minC11 + (static_cast<float>(C11[aId]) - 1.0f) * (maxC11 - minC11) / (COMPRESS - 1))*1e9f;
    float c12 = (minC12 + (static_cast<float>(C12[aId]) - 1.0f) * (maxC12 - minC12) / (COMPRESS - 1))*1e9f;
    float c13 = (minC13 + (static_cast<float>(C13[aId]) - 1.0f) * (maxC13 - minC13) / (COMPRESS - 1))*1e9f;
    float c14 = (minC14 + (static_cast<float>(C14[aId]) - 1.0f) * (maxC14 - minC14) / (COMPRESS - 1))*1e9f;
    float c15 = (minC15 + (static_cast<float>(C15[aId]) - 1.0f) * (maxC15 - minC15) / (COMPRESS - 1))*1e9f;
    float c16 = (minC16 + (static_cast<float>(C16[aId]) - 1.0f) * (maxC16 - minC16) / (COMPRESS - 1))*1e9f;

    float c22 = (minC22 + (static_cast<float>(C22[aId]) - 1.0f) * (maxC22 - minC22) / (COMPRESS - 1))*1e9f;
    float c23 = (minC23 + (static_cast<float>(C23[aId]) - 1.0f) * (maxC23 - minC23) / (COMPRESS - 1))*1e9f;
    float c24 = (minC24 + (static_cast<float>(C24[aId]) - 1.0f) * (maxC24 - minC24) / (COMPRESS - 1))*1e9f;
    float c25 = (minC25 + (static_cast<float>(C25[aId]) - 1.0f) * (maxC25 - minC25) / (COMPRESS - 1))*1e9f;
    float c26 = (minC26 + (static_cast<float>(C26[aId]) - 1.0f) * (maxC26 - minC26) / (COMPRESS - 1))*1e9f;

    float c33 = (minC33 + (static_cast<float>(C33[aId]) - 1.0f) * (maxC33 - minC33) / (COMPRESS - 1))*1e9f;
    float c34 = (minC34 + (static_cast<float>(C34[aId]) - 1.0f) * (maxC34 - minC34) / (COMPRESS - 1))*1e9f;
    float c35 = (minC35 + (static_cast<float>(C35[aId]) - 1.0f) * (maxC35 - minC35) / (COMPRESS - 1))*1e9f;
    float c36 = (minC36 + (static_cast<float>(C36[aId]) - 1.0f) * (maxC36 - minC36) / (COMPRESS - 1))*1e9f;

    float c44 = (minC44 + (static_cast<float>(C44[aId]) - 1.0f) * (maxC44 - minC44) / (COMPRESS - 1))*1e9f;
    float c45 = (minC45 + (static_cast<float>(C45[aId]) - 1.0f) * (maxC45 - minC45) / (COMPRESS - 1))*1e9f;
    float c46 = (minC46 + (static_cast<float>(C46[aId]) - 1.0f) * (maxC46 - minC46) / (COMPRESS - 1))*1e9f;

    float c55 = (minC55 + (static_cast<float>(C55[aId]) - 1.0f) * (maxC55 - minC55) / (COMPRESS - 1))*1e9f;
    float c56 = (minC56 + (static_cast<float>(C56[aId]) - 1.0f) * (maxC56 - minC56) / (COMPRESS - 1))*1e9f;

    float c66 = (minC66 + (static_cast<float>(C66[aId]) - 1.0f) * (maxC66 - minC66) / (COMPRESS - 1))*1e9f;

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
