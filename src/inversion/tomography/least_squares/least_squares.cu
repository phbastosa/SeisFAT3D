# include "least_squares.cuh"

void Least_Squares::set_parameters()
{
    set_general_parameters();
    
    set_forward_modeling();

    set_main_components();

    dx_tomo = std::stof(catch_parameter("dx_tomo", file));
    dy_tomo = std::stof(catch_parameter("dy_tomo", file));
    dz_tomo = std::stof(catch_parameter("dz_tomo", file));

    nz_tomo = (int)((modeling->nz-1) * modeling->dz / dz_tomo) + 1;    
    nx_tomo = (int)((modeling->nx-1) * modeling->dx / dx_tomo) + 1;    
    ny_tomo = (int)((modeling->ny-1) * modeling->dy / dy_tomo) + 1;  

    n_model = nx_tomo * ny_tomo * nz_tomo;
}

void Least_Squares::forward_modeling()
{
    for (int shot = 0; shot < modeling->total_shots; shot++)
    {
        modeling->shot_id = shot;

        tomography_message();

        modeling->initial_setup();
        modeling->forward_solver();
        modeling->build_outputs();

        extract_calculated_data();

        gradient_ray_tracing();
    }
}

void Least_Squares::gradient_ray_tracing()
{
    int nxx = modeling->nxx;
    int nzz = modeling->nzz;

    float rayStep = 0.2f * modeling->dz;

    for (int ray_id = 0; ray_id < modeling->total_nodes; ray_id++)
    {
        float zi = modeling->geometry->nodes.z[ray_id];
        float xi = modeling->geometry->nodes.x[ray_id];
        float yi = modeling->geometry->nodes.y[ray_id];

        std::vector < int > ray_index;

        while (true)
        {
            int i = (int)(zi / modeling->dz);
            int j = (int)(xi / modeling->dx);
            int k = (int)(yi / modeling->dy);

            float dTz = (modeling->T[(i+1) + j*nzz + k*nxx*nzz] - modeling->T[(i-1) + j*nzz + k*nxx*nzz]) / (2.0f*modeling->dz);    
            float dTx = (modeling->T[i + (j+1)*nzz + k*nxx*nzz] - modeling->T[i + (j-1)*nzz + k*nxx*nzz]) / (2.0f*modeling->dx);    
            float dTy = (modeling->T[i + j*nzz + (k+1)*nxx*nzz] - modeling->T[i + j*nzz + (k-1)*nxx*nzz]) / (2.0f*modeling->dy);

            float norm = sqrtf(dTx*dTx + dTy*dTy + dTz*dTz);

            zi -= rayStep*dTz / norm;    
            xi -= rayStep*dTx / norm;   
            yi -= rayStep*dTy / norm;   

            int im = (int)(zi / dz_tomo); 
            int jm = (int)(xi / dx_tomo); 
            int km = (int)(yi / dy_tomo); 

            ray_index.push_back(im + jm*nz_tomo + km*nx_tomo*nz_tomo);

            if (ray_index.back() == modeling->source_id) break;
        }
    
        float final_distance = sqrtf(powf(zi - modeling->geometry->shots.z[modeling->shot_id],2.0f) + 
                                     powf(xi - modeling->geometry->shots.x[modeling->shot_id],2.0f) + 
                                     powf(yi - modeling->geometry->shots.y[modeling->shot_id],2.0f));

        std::sort(ray_index.begin(), ray_index.end());

        int current = ray_index[0];
        float distance = rayStep;

        for (int index = 0; index < ray_index.size(); index++)
        {
            if (ray_index[index] == current)
            {
                distance += rayStep;
            }
            else
            {
                vG.push_back(distance);
                jG.push_back(current);
                iG.push_back(ray_id + modeling->shot_id * modeling->total_nodes);

                if (current == modeling->source_id) 
                    vG.back() = final_distance;

                distance = rayStep;
                current = ray_index[index];    
            }
        }

        if (current == modeling->source_id)
        {
            vG.push_back(final_distance);
            jG.push_back(modeling->source_id);
            iG.push_back(ray_id + modeling->shot_id * modeling->total_nodes);
        }
        else 
        {
            vG.push_back(distance);
            jG.push_back(current);
            iG.push_back(ray_id + modeling->shot_id * modeling->total_nodes);
        }

        std::vector < int >().swap(ray_index);
    }
}

void Least_Squares::optimization()
{
    std::cout<<"Solving linear system using Tikhonov regularization with order "+ std::to_string(tk_order) + "\n\n";

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
    slowness_variation_rescaling();

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
        vA[index] = lambda * vL[index - (NNZ - nnz)];        
    }

    delete[] iL;
    delete[] jL;
    delete[] vL;
}

void Least_Squares::solve_linear_system_lscg()
{
    int NIT = 10;
    float TOL = 1e-6f;

    cublasHandle_t cublasHandle = 0;
    cublasStatus_t cublasStatus;

    cublasStatus = cublasCreate(&cublasHandle);
    if (cublasStatus != CUBLAS_STATUS_SUCCESS) printf("Error cublas\n");

    cusparseHandle_t cusparseHandle = 0;
    cusparseStatus_t cusparseStatus;

    cusparseStatus = cusparseCreate(&cusparseHandle);
    if (cusparseStatus != CUSPARSE_STATUS_SUCCESS) printf("Error cusparse\n");

    size_t bsize;
    void * buffer;
    float beta = 0.0f;
    float alpha = 1.0f;
    float a, b, qTq, rTr, rd;

    int * d_iA_coo; cudaMalloc((void **)&d_iA_coo, NNZ * sizeof(int));
    int * d_iA_csr; cudaMalloc((void **)&d_iA_csr,(N+1)* sizeof(int));

    cudaMemcpy(d_iA_coo, iA, NNZ * sizeof(int), cudaMemcpyHostToDevice);

    cusparseXcoo2csr(cusparseHandle, d_iA_coo, NNZ, N, d_iA_csr, CUSPARSE_INDEX_BASE_ZERO);

    cudaFree(d_iA_coo);

    float * d_p; cudaMalloc((void **)&d_p, M * sizeof(float)); 
    float * d_q; cudaMalloc((void **)&d_q, N * sizeof(float));  
    float * d_r; cudaMalloc((void **)&d_r, M * sizeof(float)); 
    float * d_s; cudaMalloc((void **)&d_s, N * sizeof(float)); 
    float * d_x; cudaMalloc((void **)&d_x, M * sizeof(float)); 

    float * d_vA; cudaMalloc((void **)&d_vA, NNZ * sizeof(float));     
    int * d_jA_coo; cudaMalloc((void **)&d_jA_coo, NNZ * sizeof(int)); 

    cudaMemset(d_x, 0, M * sizeof(float));    
    cudaMemset(d_p, 0, M * sizeof(float));
    cudaMemset(d_q, 0, N * sizeof(float));
    cudaMemset(d_r, 0, M * sizeof(float));
    cudaMemcpy(d_s, B, N * sizeof(float), cudaMemcpyHostToDevice);

    cudaMemcpy(d_vA, vA, NNZ * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_jA_coo, jA, NNZ * sizeof(int), cudaMemcpyHostToDevice);

    cusparseDnVecDescr_t Dn_p;    
    cusparseDnVecDescr_t Dn_q;    
    cusparseDnVecDescr_t Dn_r;    
    cusparseDnVecDescr_t Dn_s;    
    cusparseSpMatDescr_t Sp_matA; 
    
    cusparseCreateDnVec(&Dn_p, M, d_p, CUDA_R_32F);
    cusparseCreateDnVec(&Dn_q, N, d_q, CUDA_R_32F);
    cusparseCreateDnVec(&Dn_r, M, d_r, CUDA_R_32F);
    cusparseCreateDnVec(&Dn_s, N, d_s, CUDA_R_32F);

    cusparseCreateCsr(&Sp_matA, N, M, NNZ, d_iA_csr, d_jA_coo, d_vA, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);

    alpha = 1.0f;
    beta = 0.0f;

    cusparseSpMV_bufferSize(cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE, &alpha, Sp_matA, Dn_s, &beta, Dn_r, CUDA_R_32F, CUSPARSE_SPMV_CSR_ALG1, &bsize);
    cudaMalloc(&buffer, bsize);

    cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE, &alpha, Sp_matA, Dn_s, &beta, Dn_r, CUDA_R_32F, CUSPARSE_SPMV_CSR_ALG1, buffer);
    cudaDeviceSynchronize();    

    cublasScopy_v2(cublasHandle, M, d_r, 1, d_p, 1);

    cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, Sp_matA, Dn_p, &beta, Dn_q, CUDA_R_32F, CUSPARSE_SPMV_CSR_ALG1, buffer);
    cudaDeviceSynchronize();    

    for (int iteration = 0; iteration < NIT; iteration++)
    {
        qTq = 0.0f;
        cublasSdot_v2(cublasHandle, N, d_q, 1, d_q, 1, &qTq);
        cudaDeviceSynchronize();    // qTq = q' * q

        rTr = 0.0f;
        cublasSdot_v2(cublasHandle, M, d_r, 1, d_r, 1, &rTr);
        cudaDeviceSynchronize();    // rTr = r' * r 

        a = rTr / qTq;              // a = (r' * r) / (q' * q)
        cublasSaxpy_v2(cublasHandle, M, &a, d_p, 1, d_x, 1);
        cudaDeviceSynchronize();    // x = x + a * p

        a *= -1.0f;
        cublasSaxpy_v2(cublasHandle, N, &a, d_q, 1, d_s, 1);
        cudaDeviceSynchronize();    // s = s - a * q 

        rd = 0.0f;
        cublasSdot_v2(cublasHandle, M, d_r, 1, d_r, 1, &rd);
        cudaDeviceSynchronize();    // rd = r' * r

        if (sqrtf(rd) < TOL) break; // Convergence condition 

        cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE, &alpha, Sp_matA, Dn_s, &beta, Dn_r, CUDA_R_32F, CUSPARSE_SPMV_CSR_ALG1, buffer);
        cudaDeviceSynchronize();   // r = G' * s    

        rTr = 0.0f;
        cublasSdot_v2(cublasHandle, M, d_r, 1, d_r, 1, &rTr);
        cudaDeviceSynchronize();   // rTr = r' * r 

        b = rTr / rd;              // b = (r' * r) / rd  
        cublasSscal_v2(cublasHandle, M, &b, d_p, 1);
        cudaDeviceSynchronize();   // p = b * p  

        b = 1.0f;
        cublasSaxpy_v2(cublasHandle, M, &b, d_r, 1, d_p, 1);
        cudaDeviceSynchronize();   // p += r  <---> p = r + b * p  

        cusparseSpMV(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, Sp_matA, Dn_p, &beta, Dn_q, CUDA_R_32F, CUSPARSE_SPMV_CSR_ALG1, buffer);
        cudaDeviceSynchronize();   // q = G * p    
    }

    cudaMemcpy(x, d_x, M * sizeof(float), cudaMemcpyDeviceToHost);

    cusparseDestroyDnVec(Dn_p);
    cusparseDestroyDnVec(Dn_q);
    cusparseDestroyDnVec(Dn_r);
    cusparseDestroyDnVec(Dn_s);
    cusparseDestroySpMat(Sp_matA);

    cudaFree(d_vA);
    cudaFree(d_iA_csr);
    cudaFree(d_jA_coo);

    cudaFree(d_x);
    cudaFree(d_p);
    cudaFree(d_q);
    cudaFree(d_r);
    cudaFree(d_s);

    cusparseDestroy(cusparseHandle);
    cublasDestroy_v2(cublasHandle);
}

void Least_Squares::slowness_variation_rescaling()
{
    for (int index = 0; index < modeling->nPoints; index++)
    {
        int k = (int) (index / (modeling->nx*modeling->nz));        
        int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;    
        int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz);  

        float xp = j*modeling->dx; 
        float yp = k*modeling->dy; 
        float zp = i*modeling->dz; 

        float x0 = floorf(xp/dx_tomo)*dx_tomo;
        float y0 = floorf(yp/dy_tomo)*dy_tomo;
        float z0 = floorf(zp/dz_tomo)*dz_tomo;

        float x1 = floorf(xp/dx_tomo)*dx_tomo + dx_tomo;
        float y1 = floorf(yp/dy_tomo)*dy_tomo + dy_tomo;
        float z1 = floorf(zp/dz_tomo)*dz_tomo + dz_tomo;

        dm[index] = 0.0f;

        if ((i >= 0) && (i < modeling->nz) && (j >= 0) && (j < modeling->nx) && (k >= 0) && (k < modeling->ny))
        {
            int idz = (int)(zp/dz_tomo);
            int idx = (int)(xp/dx_tomo);
            int idy = (int)(yp/dy_tomo);

            int ind_m = (int)(idz + idx*nz_tomo + idy*nx_tomo*nz_tomo);

            float c000 = x[ind_m];                  
            float c001 = x[ind_m + 1];
            float c100 = x[ind_m + nz_tomo];
            float c101 = x[ind_m + 1 + nz_tomo];
            float c010 = x[ind_m + nx_tomo*nz_tomo];
            float c011 = x[ind_m + 1 + nx_tomo*nz_tomo];
            float c110 = x[ind_m + nz_tomo + nx_tomo*nz_tomo];
            float c111 = x[ind_m + 1 + nz_tomo + nx_tomo*nz_tomo];  

            float xd = (xp - x0) / (x1 - x0);
            float yd = (yp - y0) / (y1 - y0);
            float zd = (zp - z0) / (z1 - z0);

            float c00 = c000*(1 - xd) + c100*xd;    
            float c01 = c001*(1 - xd) + c101*xd;    
            float c10 = c010*(1 - xd) + c110*xd;    
            float c11 = c011*(1 - xd) + c111*xd;    

            float c0 = c00*(1 - yd) + c10*yd;
            float c1 = c01*(1 - yd) + c11*yd;

            float dm_ijk = (c0*(1 - zd) + c1*zd);

            dm[i + j*modeling->nz + k*modeling->nx*modeling->nz] = dm_ijk;            
        }
    }
}