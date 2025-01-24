# include "least_squares.hpp"

void Least_Squares::set_specifications()
{
    inversion_name = "least_squares_";
    inversion_method = "Least-Squares First-Arrival Tomography";

    tk_order = std::stoi(catch_parameter("tk_order", parameters));
    tk_param = std::stof(catch_parameter("tk_param", parameters));

    n_model = modeling->nPoints;

    ray_path_max_samples = 0;

    for (int shot = 0; shot < modeling->geometry->nrel; shot++)
    {
        for (int node = modeling->geometry->iRec[shot]; node < modeling->geometry->fRec[shot]; node++)
        {
            float dx = (modeling->geometry->xsrc[modeling->geometry->sInd[shot]] - modeling->geometry->xrec[node]) / modeling->dx;
            float dy = (modeling->geometry->ysrc[modeling->geometry->sInd[shot]] - modeling->geometry->yrec[node]) / modeling->dy;
            float dz = (modeling->geometry->zsrc[modeling->geometry->sInd[shot]] - modeling->geometry->zrec[node]) / modeling->dz;
            
            ray_path_max_samples += (size_t)(sqrtf(dx*dx + dy*dy + dz*dz));
        }
    }

    iG.reserve(ray_path_max_samples);
    jG.reserve(ray_path_max_samples);
    vG.reserve(ray_path_max_samples);
}

void Least_Squares::apply_inversion_technique()
{
    float rayStep = 0.2f * modeling->dz;

    int sId = (modeling->sIdz - modeling->nb) + 
              (modeling->sIdx - modeling->nb)*modeling->nz + 
              (modeling->sIdy - modeling->nb)*modeling->nx*modeling->nz; 

    std::vector < int > ray_index; 

    for (int ray_id = modeling->geometry->iRec[modeling->srcId]; ray_id < modeling->geometry->fRec[modeling->srcId]; ray_id++)
    {
        float xi = modeling->geometry->xrec[ray_id];        
        float yi = modeling->geometry->yrec[ray_id];        
        float zi = modeling->geometry->zrec[ray_id];

        if ((modeling->geometry->xsrc[modeling->geometry->sInd[modeling->srcId]] == xi) && 
            (modeling->geometry->ysrc[modeling->geometry->sInd[modeling->srcId]] == yi) &&
            (modeling->geometry->zsrc[modeling->geometry->sInd[modeling->srcId]] == zi))
            continue;        

        while (true)
        {
            int i = (int)(zi / modeling->dz) + modeling->nb;
            int j = (int)(xi / modeling->dx) + modeling->nb;
            int k = (int)(yi / modeling->dy) + modeling->nb;

            float dTz = (modeling->T[(i+1) + j*modeling->nzz + k*modeling->nxx*modeling->nzz] - modeling->T[(i-1) + j*modeling->nzz + k*modeling->nxx*modeling->nzz]) / (2.0f*modeling->dz);    
            float dTx = (modeling->T[i + (j+1)*modeling->nzz + k*modeling->nxx*modeling->nzz] - modeling->T[i + (j-1)*modeling->nzz + k*modeling->nxx*modeling->nzz]) / (2.0f*modeling->dx);    
            float dTy = (modeling->T[i + j*modeling->nzz + (k+1)*modeling->nxx*modeling->nzz] - modeling->T[i + j*modeling->nzz + (k-1)*modeling->nxx*modeling->nzz]) / (2.0f*modeling->dy);    

            float norm = sqrtf(dTx*dTx + dTy*dTy + dTz*dTz);

            xi -= rayStep*dTx / norm;   
            yi -= rayStep*dTy / norm;   
            zi -= rayStep*dTz / norm;    

            int km = (int)(yi / modeling->dy); 
            int jm = (int)(xi / modeling->dx); 
            int im = (int)(zi / modeling->dz); 

            int index = im + jm*modeling->nz + km*modeling->nx*modeling->nz;
            
            ray_index.push_back(index);

            if (ray_index.back() == sId) break;
        }
   
        float final_distance = sqrtf(powf(xi - modeling->geometry->xsrc[modeling->geometry->sInd[modeling->srcId]], 2.0f) + 
                                     powf(yi - modeling->geometry->ysrc[modeling->geometry->sInd[modeling->srcId]], 2.0f) +
                                     powf(zi - modeling->geometry->zsrc[modeling->geometry->sInd[modeling->srcId]], 2.0f));

        std::sort(ray_index.begin(), ray_index.end());

        int current_voxel_index = ray_index[0];
        float distance_per_voxel = rayStep;

        for (int index = 0; index < ray_index.size(); index++)
        {
            if (ray_index[index] == current_voxel_index)
            {
                distance_per_voxel += rayStep;
            }
            else
            {
                vG.emplace_back(distance_per_voxel);
                jG.emplace_back(current_voxel_index);
                iG.emplace_back(ray_id + modeling->srcId * modeling->geometry->spread[modeling->srcId]);

                if (current_voxel_index == sId) vG.back() = final_distance;

                distance_per_voxel = rayStep;
                current_voxel_index = ray_index[index];    
            }
        }

        if (current_voxel_index == sId)
        {
            vG.emplace_back(final_distance);
            jG.emplace_back(current_voxel_index);
            iG.emplace_back(ray_id + modeling->srcId * modeling->geometry->spread[modeling->srcId]);
        }
        else 
        {
            vG.emplace_back(distance_per_voxel);
            jG.emplace_back(current_voxel_index);
            iG.emplace_back(ray_id + modeling->srcId * modeling->geometry->spread[modeling->srcId]);
        }

        std::vector<int>().swap(ray_index);
    }
}

void Least_Squares::optimization()
{
    M = n_model;                                  
    N = n_data + n_model - tk_order;                    
    NNZ = vG.size() + (tk_order + 1) * (n_model - tk_order);

    iA = new int[NNZ]();
    jA = new int[NNZ]();
    vA = new float[NNZ]();

    B = new float[N]();
    x = new float[M]();

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

    apply_regularization();
    solve_linear_system_lscg();

    # pragma omp parallel for
    for (int index = 0; index < n_model; index++)
        perturbation[index] = x[index];

    delete[] x;
    delete[] B;
    delete[] iA;
    delete[] jA;
    delete[] vA;
}

void Least_Squares::apply_regularization()
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

void Least_Squares::solve_linear_system_lscg()
{
    float a, b, qTq, rTr, rd;
    int cg_max_iteration = 10;

    float * s = new float[N]();
    float * q = new float[N]();
    float * r = new float[M]();
    float * p = new float[M]();

    // s = d - G * x, where d = dobs - dcal and x = slowness variation
    for (int i = 0; i < N; i++) 
        s[i] = B[i]; 

    // r = G' * s    
    for (int i = 0; i < NNZ; i++) 
        r[jA[i]] += vA[i] * s[iA[i]];        

    // p = r and x = 0;
    for (int i = 0; i < M; i++) 
    {
        p[i] = r[i]; 
        x[i] = 0.0f;
    }

    // q = G * p
    for (int i = 0; i < NNZ; i++) 
        q[iA[i]] += vA[i] * p[jA[i]];        

    for (int i = 0; i < cg_max_iteration; i++)
    {
        qTq = 0.0f;
        for (int k = 0; k < N; k++)           // q inner product
            qTq += q[k] * q[k];               // qTq = q' * q

        rTr = 0.0f;
        for (int k = 0; k < M; k++)           // r inner product
            rTr += r[k] * r[k];               // rTr = r' * r 

        a = rTr / qTq;                        // a = (r' * r) / (q' * q)                    

        for (int k = 0; k < M; k++)           // model atualization
            x[k] += a * p[k];                 // x = x + a * p

        for (int k = 0; k < N; k++)           // s atualization  
            s[k] -= a * q[k];                 // s = s - a * q 

        rd = 0.0f;
        for (int k = 0; k < M; k++)           // r inner product for division 
            rd += r[k] * r[k];                // rd = r' * r

        for (int k = 0; k < M; k++)           // Zeroing r 
            r[k] = 0.0f;                      // r = 0, for multiplication
        
        for (int k = 0; k < NNZ; k++)         // r atualization 
            r[jA[k]] += vA[k] * s[iA[k]];     // r = G' * s    

        rTr = 0.0f;                
        for (int k = 0; k < M; k++)           // r inner product
            rTr += r[k] * r[k];               // rTr = r' * r

        b = rTr / rd;                         // b = (r' * r) / rd

        for (int k = 0; k < M; k++)          
            p[k] = r[k] + b * p[k];           // p = r + b * p 

        for (int k = 0; k < N; k++) 
            q[k] = 0.0f;                      // q = 0, for multiplication

        for (int k = 0; k < NNZ; k++) 
            q[iA[k]] += vA[k] * p[jA[k]];     // q = G * p   
    }
}
