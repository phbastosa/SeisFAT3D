# include "tomography_iso.hpp"

void Tomography_ISO::set_modeling_type()
{
    modeling = new Eikonal_ISO();
    modeling->parameters = parameters;
    modeling->set_parameters();

    inversion_name = "tomography_iso";
    inversion_method = "Isotropic First-Arrival Tomography";

    dS = new float[modeling->nPoints]();
}

void Tomography_ISO::set_sensitivity_matrix()
{
    int nnz = (tk_order + 1) * (n_model - tk_order);                    

    int gsize = vG.size();

    M = n_model;                                  
    N = n_data + (n_model - tk_order);                    
    NNZ = gsize + nnz;

    iA = new int[NNZ]();
    jA = new int[NNZ]();
    vA = new float[NNZ]();

    B = new float[N]();
    x = new float[M]();
    
    for (int index = 0; index < n_data; index++)
        W[index] = 0.0f;    

    for (int index = 0; index < n_model; index++)
        R[index] = 0.0f;    

    for (int index = 0; index < gsize; index++)
    {
        W[iG[index]] += vG[index];
        R[jG[index]] += vG[index];
    }   

    for (int index = 0; index < n_data; index++) 
        B[index] = (dobs[index] - dcal[index]) * powf(1.0f/W[index], 2.0f);

    for (int index = 0; index < gsize; index++)
    {
        iA[index] = iG[index];
        jA[index] = jG[index];
        vA[index] = vG[index] * powf(1.0f/W[iG[index]], 2.0f);
    }

    for (int index = 0; index < nnz; index++)
    {
        iA[index] = iR[index] + n_data;
        jA[index] = jR[index];
        vA[index] = vR[index] * R[jR[index]]*tk_param*tk_param;  
    }

    std::vector< int >().swap(iG);
    std::vector< int >().swap(jG);
    std::vector<float>().swap(vG);        

    delete[] iR;
    delete[] jR;
    delete[] vR;
}

void Tomography_ISO::get_parameter_variation()
{
    #pragma omp parallel for
    for (int index = 0; index < modeling->nPoints; index++)
        dS[index] = x[index];
}

void Tomography_ISO::model_update()
{
    model_smoothing(dS);

    # pragma omp parallel for
    for (int index = 0; index < modeling->nPoints; index++)
    {
        int k = (int) (index / (modeling->nx*modeling->nz));         
        int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;   
        int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz); 

        int indb = (i + modeling->nb) + (j + modeling->nb)*modeling->nzz;

        # pragma omp atomic
        modeling->S[indb] += dS[index];
    }

    modeling->copy_slowness_to_device();
}

void Tomography_ISO::export_estimated_models()
{
    float * Vp = new float[modeling->nPoints]();
    modeling->reduce_boundary(modeling->S, Vp);
    
    # pragma omp parallel for
    for (int index = 0; index < modeling->nPoints; index++)
        Vp[index] = 1.0f / Vp[index];

    std::string estimated_vp_path = estimated_model_folder + inversion_name + "_final_model_vp_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + ".bin";
    export_binary_float(estimated_vp_path, Vp, modeling->nPoints);
}
