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

    std::string Cijkl_folder = catch_parameter("Cijkl_folder", parameters);

    C11 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C11.bin", C11, nPoints);

    C12 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C12.bin", C12, nPoints);

    C13 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C13.bin", C13, nPoints);

    C14 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C14.bin", C14, nPoints);

    C15 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C15.bin", C15, nPoints);

    C16 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C16.bin", C16, nPoints);

    C22 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C22.bin", C22, nPoints);

    C23 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C23.bin", C23, nPoints);
    
    C24 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C24.bin", C24, nPoints);

    C25 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C25.bin", C25, nPoints);

    C26 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C26.bin", C26, nPoints);

    C33 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C33.bin", C33, nPoints);
    
    C34 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C34.bin", C34, nPoints);

    C35 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C35.bin", C35, nPoints);

    C36 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C36.bin", C36, nPoints);

    C44 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C44.bin", C44, nPoints);

    C45 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C45.bin", C45, nPoints);

    C46 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C46.bin", C46, nPoints);

    C55 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C55.bin", C55, nPoints);

    C56 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C56.bin", C56, nPoints);

    C66 = new float[nPoints]();
    import_binary_float(Cijkl_folder + "C66.bin", C66, nPoints);
}

void Eikonal_ANI::forward_solver()
{
    cudaMemcpy(d_S, S, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_T, T, volsize*sizeof(float), cudaMemcpyHostToDevice);

    fast_sweeping_method();

    cudaMemcpy(T, d_T, volsize*sizeof(float), cudaMemcpyDeviceToHost);

    int sIdx = (int)(geometry->xrec[geometry->sInd[srcId]] / dx) + nb;
    int sIdy = (int)(geometry->xrec[geometry->sInd[srcId]] / dy) + nb;
    int sIdz = (int)(geometry->xrec[geometry->sInd[srcId]] / dz) + nb;

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

            p[0] = dTz / norm;
            p[1] = dTx / norm;
            p[2] = dTy / norm;

            get_stiffness();
            get_christoffel();
            get_eigen_values();

            qS[index] = 1.0f / sqrtf(Gv[0] * Ro[aId]);
        }
    }

    initialization();

    cudaMemcpy(d_S, qS, volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_T, T, volsize*sizeof(float), cudaMemcpyHostToDevice);

    fast_sweeping_method();

    cudaMemcpy(T, d_T, volsize*sizeof(float), cudaMemcpyDeviceToHost);

    compute_seismogram();
}

void Eikonal_ANI::get_stiffness()
{
    C[0+0*v] = C11[aId]; C[0+1*v] = C12[aId]; C[0+2*v] = C13[aId]; 
    C[1+0*v] = C12[aId]; C[1+1*v] = C22[aId]; C[1+2*v] = C23[aId]; 
    C[2+0*v] = C13[aId]; C[2+1*v] = C23[aId]; C[2+2*v] = C33[aId]; 
    C[3+0*v] = C14[aId]; C[3+1*v] = C24[aId]; C[3+2*v] = C34[aId]; 
    C[4+0*v] = C15[aId]; C[4+1*v] = C25[aId]; C[4+2*v] = C35[aId]; 
    C[5+0*v] = C16[aId]; C[5+1*v] = C26[aId]; C[5+2*v] = C36[aId]; 
    
    C[0+3*v] = C14[aId]; C[0+4*v] = C15[aId]; C[0+5*v] = C16[aId];
    C[1+3*v] = C24[aId]; C[1+4*v] = C25[aId]; C[1+5*v] = C26[aId];
    C[2+3*v] = C34[aId]; C[2+4*v] = C35[aId]; C[2+5*v] = C36[aId];
    C[3+3*v] = C44[aId]; C[3+4*v] = C45[aId]; C[3+5*v] = C46[aId];
    C[4+3*v] = C45[aId]; C[4+4*v] = C55[aId]; C[4+5*v] = C56[aId];
    C[5+3*v] = C46[aId]; C[5+4*v] = C56[aId]; C[5+5*v] = C66[aId];

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
    if (Gv[0] > Gv[1]) std::swap(Gv[0],Gv[1]);    
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
