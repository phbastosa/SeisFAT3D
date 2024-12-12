# include "eikonal_vti.cuh"

void Eikonal_VTI::set_properties()
{
    S = new float[volsize]();

    Vp = new float[nPoints]();
    Vs = new float[nPoints]();
    Ro = new float[nPoints]();

    E = new float[nPoints]();
    D = new float[nPoints]();
    G = new float[nPoints](); 

    S_vti = new float[volsize]();

    std::string vp_file = catch_parameter("vp_model_file", parameters);
    std::string vs_file = catch_parameter("vs_model_file", parameters);
    std::string ro_file = catch_parameter("ro_model_file", parameters);

    std::string e_file = catch_parameter("epsilon_model_file", parameters);
    std::string d_file = catch_parameter("delta_model_file", parameters);
    std::string g_file = catch_parameter("gamma_model_file", parameters);

    import_binary_float(vp_file, Vp, nPoints);
    import_binary_float(vs_file, Vs, nPoints);
    import_binary_float(ro_file, Ro, nPoints);

    import_binary_float(e_file, E, nPoints);
    import_binary_float(d_file, D, nPoints);
    import_binary_float(g_file, G, nPoints);
    
    float * slowness = new float[nPoints]();

    for (int index = 0; index < nPoints; index++)
        slowness[index] = 1.0f / Vp[index];

    expand_boundary(slowness, S);

    for (int index = 0; index < volsize; index++)
        S_vti[index] = S[index];

    delete[] slowness;
}

void Eikonal_VTI::set_conditions()
{
    modeling_type = "eikonal_vti";
    modeling_name = "Modeling type: Eikonal vertically transverse isotropic time propagation";

    n = 3;
    v = 6;

    p = new float[n]();
    Gv = new float[n]();
    Gij = new float[n*n]();
    Cijkl = new float[v*v]();
}

void Eikonal_VTI::forward_solver()
{
    cudaMemcpy(d_S, S, volsize*sizeof(float), cudaMemcpyHostToDevice);

    initialization();

    cudaMemcpy(d_T, T, volsize*sizeof(float), cudaMemcpyHostToDevice);

    fast_sweeping_method();

    cudaMemcpy(T, d_T, volsize*sizeof(float), cudaMemcpyDeviceToHost);

    int sIdx = (int)(geometry->xrec[geometry->sInd[srcId]] / dx);
    int sIdy = (int)(geometry->xrec[geometry->sInd[srcId]] / dy);
    int sIdz = (int)(geometry->xrec[geometry->sInd[srcId]] / dz);

    for (int index = 0; index < nPoints; index++)
    {
        int k = (int) (index / (nx*nz));         
        int j = (int) (index - k*nx*nz) / nz;    
        int i = (int) (index - j*nz - k*nx*nz);  

        if ((i == sIdz) && (j == sIdx) && (k == sIdy))    
            continue;

        aId = i + j*nx + k*nx*nz;
        bId = (i + nb) + (j + nb)*nzz + (k + nb)*nxx*nzz;

        float dTz = 0.5f*(T[bId + 1] - T[bId - 1]) / dz;
        float dTx = 0.5f*(T[bId + nzz] - T[bId - nzz]) / dx;
        float dTy = 0.5f*(T[bId + nxx*nzz] - T[bId - nxx*nzz]) / dy;

        float norm = sqrtf(dTx*dTx + dTy*dTy + dTz*dTz);

        p[0] = dTz / norm;
        p[1] = dTx / norm;
        p[2] = dTy / norm;

        get_stiffness();    
        get_christoffel();
        get_eigen_values();

        S_vti[bId] = 1.0f / sqrtf(Gv[0] * Ro[aId]);
    }

    cudaMemcpy(d_S, S_vti, volsize*sizeof(float), cudaMemcpyHostToDevice);

    initialization();

    cudaMemcpy(d_T, T, volsize*sizeof(float), cudaMemcpyHostToDevice);

    fast_sweeping_method();

    cudaMemcpy(T, d_T, volsize*sizeof(float), cudaMemcpyDeviceToHost);

    compute_seismogram();
}

int Eikonal_VTI::voigt_map(int i, int j)
{
    if (i == j)
        return i;
    if (((i == 1) && (j == 2)) || ((i == 2) && (j == 1)))
        return 3;
    if (((i == 2) && (j == 0)) || ((i == 0) && (j == 2)))
        return 4;
    if (((i == 0) && (j == 1)) || ((i == 1) && (j == 0)))
        return 5;
    
    return -9999999;
}

void Eikonal_VTI::get_stiffness()
{
    float A33 = Vp[aId]*Vp[aId]*Ro[aId];
    float A55 = Vs[aId]*Vs[aId]*Ro[aId];
    float A11 = A33*(1 + 2*E[aId]); 
    float A66 = A55*(1 + 2*G[aId]);
    float A13 = sqrtf(2*D[aId]*A33*(A33 - A55) + (A33 - A55)*(A33 - A55)) - A55;

    float B2 = 1.0f / Ro[aId] / Ro[aId];

    Cijkl[2 + 2*v] = A33 * B2;
    Cijkl[5 + 5*v] = A66 * B2;
    Cijkl[0 + 0*v] = Cijkl[1 + 1*v] = A11 * B2;
    Cijkl[0 + 2*v] = Cijkl[0 + 2*v] = A13 * B2;
    Cijkl[2 + 1*v] = Cijkl[1 + 2*v] = A13 * B2;
    Cijkl[3 + 3*v] = Cijkl[4 + 4*v] = A55 * B2;
    Cijkl[1 + 0*v] = Cijkl[0 + 1*v] = (A11 - 2.0f*A66) * B2;
}

void Eikonal_VTI::get_christoffel()
{
    for (int index = 0; index < n*n; index++) 
        Gij[index] = 0.0f; 

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

                    Gij[i + j*n] += Cijkl[I + J*v]*p[k]*p[l];
                }
            }
        }
    }
}

void Eikonal_VTI::get_eigen_values()
{
    float a = -(Gij[0] + Gij[4] + Gij[8]);
    
    float b = Gij[0]*Gij[4] + Gij[4]*Gij[8] + 
              Gij[0]*Gij[8] - Gij[3]*Gij[1] - 
              Gij[6]*Gij[6] - Gij[7]*Gij[5];
    
    float c = -(Gij[0]*(Gij[4]*Gij[8] - Gij[7]*Gij[5]) -
                Gij[3]*(Gij[1]*Gij[8] - Gij[7]*Gij[6]) +
                Gij[6]*(Gij[1]*Gij[5] - Gij[4]*Gij[6]));

    float p = b - (a*a)/3.0f;
    float q = (2.0f*a*a*a)/27.0f - (a*b)/3.0f + c;

    float detM = 0.25f*(q*q) + (p*p*p)/27.0f;

    if (detM > 0) 
    {
        float u = cbrtf(-0.5f*q + sqrtf(detM));
        float v = cbrtf(-0.5f*q - sqrtf(detM));
        
        Gv[0] = u + v - a/3.0f;
    } 
    else if (detM == 0) 
    {       
        float u = cbrt(-0.5f*q);

        Gv[0] = 2.0f*u - a/3.0f;
        Gv[1] = -u - a/3.0f;         
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
    
    if (G[0] < G[1]) std::swap(G[0],G[1]);
    if (G[1] < G[2]) std::swap(G[1],G[2]);
    if (G[0] > G[1]) std::swap(G[0],G[1]);    
}
