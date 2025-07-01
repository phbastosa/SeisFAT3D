# include "tomography_iso.hpp"

void Tomography_ISO::set_modeling_type()
{
    modeling = new Eikonal_ISO();
    modeling->parameters = parameters;
    modeling->set_parameters();

    inversion_name = "tomography_iso";
    inversion_method = "Isotropic First-Arrival Tomography";
}

void Tomography_ISO::set_objective_function()
{
    float square_difference = 0.0f;
    
    for (int i = 0; i < n_data; i++)
        square_difference += powf(dobs[i] - dcal[i], 2.0f);
    
    residuo.push_back(sqrtf(square_difference));
}

void Tomography_ISO::set_sensitivity_matrix()
{
    M = n_model;                                  
    N = n_data + n_model - tk_order;                    
    NNZ = vG.size() + (tk_order + 1) * (M - tk_order);

    iA = new int[NNZ]();
    jA = new int[NNZ]();
    vA = new float[NNZ]();

    B = new float[N]();

    for (int index = 0; index < n_data; index++) 
        B[index] = dobs[index] - dcal[index];

    for (int index = 0; index < vG.size(); index++)
    {
        iA[index] = iG[index];
        jA[index] = jG[index];
        vA[index] = vG[index];
    }

    std::vector< int >().swap(iG);
    std::vector< int >().swap(jG);
    std::vector<float>().swap(vG);        
}

void Tomography_ISO::set_regularization()
{
    int elements = tk_order + 1;
		
    int n = n_model - tk_order;
    int nnz = elements * n;	
    
    int * iL = new int[nnz]();
    int * jL = new int[nnz]();
    float * vL = new float[nnz]();

    if (tk_order <= 0)
    {
	    for (int index = 0; index < nnz; index++)
        {
            iL[index] = index;
	        jL[index] = index;
	        vL[index] = 1.0f;
	    }
    } 
    else
    {
        int * df = new int[elements]();	
        int * df1 = new int[elements + 1]();
        int * df2 = new int[elements + 1]();
        
        df[0] = -1; df[1] = 1;
        
        for (int index = 1; index < tk_order; index++)
        {
            for (int k = 0; k < elements; k++)
            {
                df2[k] = df[k];
                df1[k + 1] = df[k];

                df[k] = df1[k] - df2[k]; 
            }		 
        }
        
        for (int index = 0; index < n; index++)
        {
            for (int k = 0; k < elements; k++)
            {
                iL[elements*index + k] = index;	
                jL[elements*index + k] = index + k;
                vL[elements*index + k] = df[k];
            }	
        }

        delete[] df;
        delete[] df1;
        delete[] df2;
    }

    for (int index = NNZ - nnz; index < NNZ; index++) 
    {
        iA[index] = n_data + iL[index - (NNZ - nnz)];
        jA[index] = jL[index - (NNZ - nnz)];
        vA[index] = tk_param * vL[index - (NNZ - nnz)];        
    }

    delete[] iL;
    delete[] jL;
    delete[] vL;
}