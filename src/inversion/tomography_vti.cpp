# include "tomography_vti.hpp"

void Tomography_VTI::set_modeling_type()
{
    modeling = new Eikonal_ANI();
    modeling->parameters = parameters;
    modeling->set_parameters();

    eikonal = dynamic_cast<Eikonal_ANI*>(modeling);

    inversion_name = "tomography_vti";
    inversion_method = "Anisotropic First-Arrival Tomography";

    E = new float[modeling->nPoints]();    
    D = new float[modeling->nPoints]();    

    dS = new float[modeling->nPoints]();
    dE = new float[modeling->nPoints]();    
    dD = new float[modeling->nPoints]();    

    eikonal->get_stiffness_VTI(E,D);    
}

void Tomography_VTI::set_sensitivity_matrix()
{
    int np = 3;
    int gsize = vG.size();

    int n = (n_model - tk_order);
    int nnz = (tk_order + 1) * n;   
    
    M = np*n_model;                                  
    N = n_data + np*n;
    NNZ = np*(gsize + nnz);    

    iA = new int[NNZ]();
    jA = new int[NNZ]();
    vA = new float[NNZ]();

    B = new float[N]();
    x = new float[M]();    

    # pragma omp parallel for
    for (int index = 0; index < n_data; index++)
        W[index] = 0.0f;    

    # pragma omp parallel for
    for (int index = 0; index < n_model; index++)
        R[index] = 0.0f;    

    # pragma omp parallel for        
    for (int index = 0; index < gsize; index++)
    {
        # pragma omp atomic
        W[iG[index]] += vG[index];
        
        # pragma omp atomic
        R[jG[index]] += vG[index];
    }   

    # pragma omp parallel for
    for (int index = 0; index < n_data; index++) 
        B[index] = (dobs[index] - dcal[index]) * powf(1.0f/W[index], 2.0f);

    # pragma omp parallel for
    for (int index = 0; index < gsize; index++)
    {
        int k = (int) (iG[index] / (modeling->nx*modeling->nz)) + modeling->nb;         
        int j = (int) (iG[index] - k*modeling->nx*modeling->nz) / modeling->nz + modeling->nb;   
        int i = (int) (iG[index] - j*modeling->nz - k*modeling->nx*modeling->nz) + modeling->nb; 

        float S = modeling->S[i + j*modeling->nzz];

        int indp = (i - modeling->nb) + (j - modeling->nb)*modeling->nz;

        float dTz = 0.5f*(modeling->T[(i+1) + j*modeling->nzz + k*modeling->nxx*modeling->nzz] - modeling->T[(i-1) + j*modeling->nzz + k*modeling->nxx*modeling->nzz]) / modeling->dz;    
        float dTx = 0.5f*(modeling->T[i + (j+1)*modeling->nzz + k*modeling->nxx*modeling->nzz] - modeling->T[i + (j-1)*modeling->nzz + k*modeling->nxx*modeling->nzz]) / modeling->dx;    
        float dTy = 0.5f*(modeling->T[i + j*modeling->nzz + (k+1)*modeling->nxx*modeling->nzz] - modeling->T[i + j*modeling->nzz + (k-1)*modeling->nxx*modeling->nzz]) / modeling->dy;    

        float norm = sqrtf(dTx*dTx + dTy*dTy + dTz*dTz);

        float px = dTx/norm;
        float py = dTy/norm;
        float pz = dTz/norm;

        float pnorm = sqrtf(px*px + py*py + pz*pz);

        float theta = acos(pz/pnorm);
        
        float sin2 = sinf(theta)*sinf(theta);
        float cos2 = cosf(theta)*cosf(theta);

        float denom = 1.0f + D[indp]*sin2*cos2 + E[indp]*sin2*sin2;

        float dqSdS = 1.0f / denom;
        float dqSdE =-(S*sin2*sin2) / (denom*denom);
        float dqSdD =-(S*sin2*cos2) / (denom*denom);

        iA[index] = iG[index];
        jA[index] = jG[index];
        vA[index] = vG[index]*dqSdS * sqrtf(1.0f / W[iG[index]]);

        iA[index + gsize] = iG[index];
        jA[index + gsize] = jG[index] + n_model;
        vA[index + gsize] = vG[index]*dqSdE * sqrtf(1.0f / W[iG[index]]);

        iA[index + 2*gsize] = iG[index];
        jA[index + 2*gsize] = jG[index] + 2*n_model;
        vA[index + 2*gsize] = vG[index]*dqSdD * sqrtf(1.0f / W[iG[index]]);
    }

    for (int index = 0; index < nnz; index++)
    {
        iA[index + np*gsize] = iR[index] + n_data;
        jA[index + np*gsize] = jR[index];
        vA[index + np*gsize] = vR[index] * R[jR[index]]*tk_param*tk_param;  
    
        iA[index + np*gsize + nnz] = iR[index] + n_data + n;
        jA[index + np*gsize + nnz] = jR[index] + n_model;
        vA[index + np*gsize + nnz] = vR[index] * R[jR[index]]*tk_param*tk_param;  

        iA[index + np*gsize + 2*nnz] = iR[index] + n_data + 2*n;
        jA[index + np*gsize + 2*nnz] = jR[index] + 2*n_model;
        vA[index + np*gsize + 2*nnz] = vR[index] * R[jR[index]]*tk_param*tk_param;      
    }
    
    std::vector< int >().swap(iG);
    std::vector< int >().swap(jG);
    std::vector<float>().swap(vG);        

    delete[] iR;
    delete[] jR;
    delete[] vR;
}

void Tomography_VTI::get_parameter_variation()
{
    #pragma omp parallel for
    for (int index = 0; index < modeling->nPoints; index++)
    {   
        dS[index] = x[index];
        dE[index] = x[index + n_model];
        dD[index] = x[index + 2*n_model];
    }
}

void Tomography_VTI::model_update()
{
    model_smoothing(dS);
    model_smoothing(dE);
    model_smoothing(dD);
    
    for (int index = 0; index < modeling->nPoints; index++)
    {
        int i = (int) (index % modeling->nz);    
        int j = (int) (index / modeling->nz); 

        int indb = (i + modeling->nb) + (j + modeling->nb)*modeling->nzz;

        modeling->S[indb] += dS[index];

        E[index] += dE[index];
        D[index] += dD[index];
    }

    modeling->copy_slowness_to_device();

    eikonal->set_stiffness_VTI(E,D);
}

void Tomography_VTI::export_estimated_models()
{
    float * V = new float[modeling->nPoints]();
    modeling->reduce_boundary(modeling->S, V);

    # pragma omp parallel for
    for (int index = 0; index < modeling->nPoints; index++)
        V[index] = 1.0f / V[index];

    std::string estimated_v_path = estimated_model_folder + inversion_name + "_final_model_vp_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + ".bin";
    std::string estimated_e_path = estimated_model_folder + inversion_name + "_final_model_ep_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + ".bin";
    std::string estimated_d_path = estimated_model_folder + inversion_name + "_final_model_dl_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + ".bin";

    export_binary_float(estimated_v_path, V, modeling->nPoints);
    export_binary_float(estimated_e_path, E, modeling->nPoints);
    export_binary_float(estimated_d_path, D, modeling->nPoints);
}