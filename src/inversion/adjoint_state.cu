# include "adjoint_state.cuh"

int Adjoint_State::iDivUp(int a, int b) 
{ 
    return ( (a % b) != 0 ) ? (a / b + 1) : (a / b); 
}

void Adjoint_State::set_specifications()
{
    inversion_name = "adjoint_state_";
    inversion_method = "Adjoint-State First-Arrival Tomography";

    adam_rate = std::stof(catch_parameter("adam_rate", parameters));

    aperture_x = std::stof(catch_parameter("inv_aperture_x", parameters));
    aperture_y = std::stof(catch_parameter("inv_aperture_y", parameters));

    nSweeps = 8;
    meshDim = 3;

    cell_area = modeling->dx*modeling->dy*modeling->dz;

    total_levels = (modeling->nxx - 1) + (modeling->nyy - 1) + (modeling->nzz - 1);

    m = new float[modeling->nPoints]();
    v = new float[modeling->nPoints]();

    h_source_grad = new float[modeling->volsize]();
    h_source_comp = new float[modeling->volsize]();
    
    h_adjoint_grad = new float[modeling->volsize]();
    h_adjoint_comp = new float[modeling->volsize]();

    gradient = new float[modeling->nPoints]();

    cudaMalloc((void**)&(d_T), modeling->volsize*sizeof(float));
    
    cudaMalloc((void**)&(d_source_grad), modeling->volsize*sizeof(float));
    cudaMalloc((void**)&(d_source_comp), modeling->volsize*sizeof(float));
    
    cudaMalloc((void**)&(d_adjoint_grad), modeling->volsize*sizeof(float));
    cudaMalloc((void**)&(d_adjoint_comp), modeling->volsize*sizeof(float));
}

void Adjoint_State::apply_inversion_technique()
{
    initialization();

    cudaMemcpy(d_T, modeling->T, modeling->volsize*sizeof(float), cudaMemcpyHostToDevice);
    
    cudaMemcpy(d_source_grad, h_source_grad, modeling->volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_source_comp, h_source_comp, modeling->volsize*sizeof(float), cudaMemcpyHostToDevice);
    
    cudaMemcpy(d_adjoint_grad, h_adjoint_grad, modeling->volsize*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_adjoint_comp, h_adjoint_comp, modeling->volsize*sizeof(float), cudaMemcpyHostToDevice);

    for (int sweep_iter = 0; sweep_iter < meshDim; sweep_iter++)
    {
        for (int sweep = 0; sweep < nSweeps; sweep++)
        { 
            int start = (sweep == 3 || sweep == 5 || sweep == 6 || sweep == 7) ? total_levels : meshDim;
            int end = (start == meshDim) ? total_levels + 1 : meshDim - 1;
            int incr = (start == meshDim) ? true : false;
        
            int xSweepOff = (sweep == 3 || sweep == 4) ? modeling->nxx : 0;
            int ySweepOff = (sweep == 2 || sweep == 5) ? modeling->nyy : 0;
            int zSweepOff = (sweep == 1 || sweep == 6) ? modeling->nzz : 0;
            
            for (int level = start; level != end; level = (incr) ? level + 1 : level - 1)
            {			
                int xs = max(1, level - (modeling->nyy + modeling->nzz));	
                int ys = max(1, level - (modeling->nxx + modeling->nzz));
            
                int xe = min(modeling->nxx, level - (meshDim - 1));
                int ye = min(modeling->nyy, level - (meshDim - 1));	
                
                int xr = xe - xs + 1;
                int yr = ye - ys + 1;
            
                int nThreads = xr * yr;
                    
                dim3 bs(16, 16, 1);
            
                if (nThreads < 256) { bs.x = xr; bs.y = yr; } 
            
                dim3 gs(iDivUp(xr, bs.x), iDivUp(yr , bs.y), 1);
            
                adjoint_state_kernel<<<gs,bs>>>(d_T, d_adjoint_grad, d_adjoint_comp, d_source_grad, d_source_comp, level, xs, ys, xSweepOff, ySweepOff, 
                                                zSweepOff, modeling->nxx, modeling->nyy, modeling->nzz, modeling->dx, modeling->dy, modeling->dz);
            
                cudaDeviceSynchronize();
            }
        }
    }

    cudaMemcpy(h_adjoint_grad, d_adjoint_grad, modeling->volsize*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_adjoint_comp, d_adjoint_comp, modeling->volsize*sizeof(float), cudaMemcpyDeviceToHost);

    float grad_max = 0.0f;
    float comp_max = 0.0f;

    for (int index = 0; index < modeling->volsize; index++)
    {
        grad_max = std::max(grad_max, h_adjoint_grad[index]);
        comp_max = std::max(comp_max, h_adjoint_comp[index]);
    }

    float cmp_x; 
    float cmp_y;

    float max_offset = 0.0f;

    int ri = modeling->geometry->iRec[modeling->srcId];
    int rf = modeling->geometry->fRec[modeling->srcId];

    float sx = modeling->geometry->xsrc[modeling->geometry->sInd[modeling->srcId]]; 
    float sy = modeling->geometry->ysrc[modeling->geometry->sInd[modeling->srcId]]; 

    for (int rId = ri; rId < rf; rId++)
    {
        float rx = modeling->geometry->xrec[rId];
        float ry = modeling->geometry->yrec[rId];

        float offset = sqrtf((sx - rx)*(sx - rx) + (sy - ry)*(sy - ry));

        if (max_offset < offset)
        {
            cmp_x = sx + 0.5f*(rx - sx);
            cmp_y = sy + 0.5f*(ry - sy);

            max_offset = offset;
        }
    }

    for (int index = 0; index < modeling->nPoints; index++) 
    {
        int k = (int) (index / (modeling->nx*modeling->nz));        
        int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;    
        int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz);  

        int indp = i + j*modeling->nz + k*modeling->nx*modeling->nz; 
        int indb = (i + modeling->nb) + (j + modeling->nb)*modeling->nzz + (k + modeling->nb)*modeling->nxx*modeling->nzz;

        h_adjoint_grad[indb] *= 1.0f / grad_max;
        h_adjoint_comp[indb] *= 1.0f / comp_max;

        float sigma_x = tanf(aperture_x * PI / 180.0f)*i*modeling->dz;
        float sigma_y = tanf(aperture_y * PI / 180.0f)*i*modeling->dz;

        float par_x = powf((j*modeling->dx - cmp_x)/(sigma_x + 1e-6f), 2.0f);
        float par_y = powf((k*modeling->dy - cmp_y)/(sigma_y + 1e-6f), 2.0f);

        float value = expf(-0.5*(par_x + par_y));

        gradient[indp] += value*(h_adjoint_grad[indb] / (h_adjoint_comp[indb] + 1e-3f)*cell_area / modeling->geometry->nrel);    
    }
}

void Adjoint_State::initialization()
{
    # pragma omp parallel for
    for (int index = 0; index < modeling->volsize; index++) 
    {
        h_source_grad[index] = 0.0f;    
        h_source_comp[index] = 0.0f;    
        
        h_adjoint_grad[index] = 1e6f;
        h_adjoint_comp[index] = 1e6f;

        int k = (int) (index / (modeling->nxx*modeling->nzz));         
        int j = (int) (index - k*modeling->nxx*modeling->nzz) / modeling->nzz;   
        int i = (int) (index - j*modeling->nzz - k*modeling->nxx*modeling->nzz); 

        if ((i == 0) || (i == modeling->nzz - 1) || 
            (j == 0) || (j == modeling->nxx - 1) || 
            (k == 0) || (k == modeling->nyy - 1))     
        {
            h_adjoint_grad[index] = 0.0f;        
            h_adjoint_comp[index] = 0.0f;        
        }    
    }

    int sId = modeling->geometry->sInd[modeling->srcId];

    int skipped = modeling->srcId * modeling->geometry->spread[modeling->srcId];

    int sIdx = (int)(modeling->geometry->xsrc[sId] / modeling->dx) + modeling->nb;
    int sIdy = (int)(modeling->geometry->ysrc[sId] / modeling->dy) + modeling->nb;
    int sIdz = (int)(modeling->geometry->zsrc[sId] / modeling->dz) + modeling->nb;

    float Sref = modeling->S[sIdz + sIdx*modeling->nzz + sIdy*modeling->nxx*modeling->nzz];    

    int spread = 0;

    for (modeling->recId = modeling->geometry->iRec[modeling->srcId]; modeling->recId < modeling->geometry->fRec[modeling->srcId]; modeling->recId++)
    {
        int rIdx = (int)(modeling->geometry->xrec[modeling->recId] / modeling->dx) + modeling->nb;
        int rIdy = (int)(modeling->geometry->yrec[modeling->recId] / modeling->dy) + modeling->nb;
        int rIdz = (int)(modeling->geometry->zrec[modeling->recId] / modeling->dz) + modeling->nb;

        for (int k = 0; k < 3; k++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int i = 0; i < 3; i++)
                {
                    int yi = rIdy + (k - 1);
                    int xi = rIdx + (j - 1);
                    int zi = rIdz + (i - 1);

                    int index = zi + xi*modeling->nzz + yi*modeling->nxx*modeling->nzz; 

                    float X = sqrtf(powf((sIdy - yi)*modeling->dy, 2.0f) + 
                                    powf((sIdx - xi)*modeling->dx, 2.0f) + 
                                    powf((sIdz - zi)*modeling->dz, 2.0f));

                    h_source_grad[index] += (dobs[spread + skipped] - modeling->T[index]) / cell_area;    
                    h_source_comp[index] += 1.0f / (X*X*Sref);
                }
            }
        }

        ++spread;
    }   
} 

void Adjoint_State::optimization()  
{
    float gdot = 0.0f;
    #pragma omp parallel for reduction(+:gdot)
    for (int index = 0; index < modeling->nPoints; index++)
        gdot += gradient[index]*gradient[index];
    
    float beta1 = 0.9f;
    float beta2 = 0.999f;

    float epsilon = 1e-8f;

    for (int index = 0; index < modeling->nPoints; index++)
    {
        gradient[index] *= 1.0f / gdot;

        m[index] = beta1*m[index] + (1.0f - beta1)*gradient[index];
        
        v[index] = beta2*v[index] + (1.0f - beta2)*gradient[index]*gradient[index];

        float m_hat = m[index] / (1.0f - powf(beta1, iteration));
        
        float v_hat = v[index] / (1.0f - powf(beta2, iteration));

        perturbation[index] = adam_rate*m_hat/(sqrtf(v_hat) + epsilon);
    }

    memset(gradient, 0.0f, modeling->nPoints);
}

__global__ void adjoint_state_kernel(float * T, float * adjoint_grad, float * adjoint_comp, float * source_grad, float * source_comp, int level, int xOffset, 
                                    int yOffset, int xSweepOffset, int ySweepOffset, int zSweepOffset, int nxx, int nyy, int nzz, float dx, float dy, float dz)
{
    int x = (blockIdx.x * blockDim.x + threadIdx.x) + xOffset;
    int y = (blockIdx.y * blockDim.y + threadIdx.y) + yOffset;

    if ((x < nxx) && (y < nyy)) 
    {
        int z = level - (x + y);
		
        if ((z >= 0) && (z < nzz))	
        {
            int i = (int)abs(z - zSweepOffset);
            int j = (int)abs(x - xSweepOffset);
            int k = (int)abs(y - ySweepOffset);

            if ((i > 0) && (i < nzz - 1) && (j > 0) && (j < nxx - 1) && (k > 0) && (k < nyy - 1))
            {
                float a1 = -1.0f*(T[i + j*nzz + k*nxx*nzz] - T[i + (j-1)*nzz + k*nxx*nzz]) / dx;                                                
                float ap1 = (a1 + fabsf(a1)) / 2.0f;
                float am1 = (a1 - fabsf(a1)) / 2.0f;

                float a2 = -1.0f*(T[i + (j+1)*nzz + k*nxx*nzz] - T[i + j*nzz + k*nxx*nzz]) / dx;
                float ap2 = (a2 + fabsf(a2)) / 2.0f;
                float am2 = (a2 - fabsf(a2)) / 2.0f;

                float b1 = -1.0f*(T[i + j*nzz + k*nxx*nzz] - T[i + j*nzz + (k-1)*nxx*nzz]) / dy;
                float bp1 = (b1 + fabsf(b1)) / 2.0f;
                float bm1 = (b1 - fabsf(b1)) / 2.0f;

                float b2 = -1.0f*(T[i + j*nzz + (k+1)*nxx*nzz] - T[i + j*nzz + k*nxx*nzz]) / dy;
                float bp2 = (b2 + fabsf(b2)) / 2.0f;
                float bm2 = (b2 - fabsf(b2)) / 2.0f;

                float c1 = -1.0f*(T[i + j*nzz + k*nxx*nzz] - T[(i-1) + j*nzz + k*nxx*nzz]) / dz;
                float cp1 = (c1 + fabsf(c1)) / 2.0f;
                float cm1 = (c1 - fabsf(c1)) / 2.0f;

                float c2 = -1.0f*(T[(i+1) + j*nzz + k*nxx*nzz] - T[i + j*nzz + k*nxx*nzz]) / dz;        
                float cp2 = (c2 + fabsf(c2)) / 2.0f;
                float cm2 = (c2 - fabsf(c2)) / 2.0f;

                float d = (ap2 - am1)/dx + (bp2 - bm1)/dy + (cp2 - cm1)/dz;

                if (fabsf(d) < 1e-6f)
                {
                    adjoint_grad[i + j*nzz + k*nxx*nzz] = 0.0f;    
                    adjoint_comp[i + j*nzz + k*nxx*nzz] = 0.0f;    
                }
                else
                {
                    float eg = (ap1*adjoint_grad[i + (j-1)*nzz + k*nxx*nzz] - am2*adjoint_grad[i + (j+1)*nzz + k*nxx*nzz]) / dx +
                               (bp1*adjoint_grad[i + j*nzz + (k-1)*nxx*nzz] - bm2*adjoint_grad[i + j*nzz + (k+1)*nxx*nzz]) / dy +
                               (cp1*adjoint_grad[(i-1) + j*nzz + k*nxx*nzz] - cm2*adjoint_grad[(i+1) + j*nzz + k*nxx*nzz]) / dz;

                    float ec = (ap1*adjoint_comp[i + (j-1)*nzz + k*nxx*nzz] - am2*adjoint_comp[i + (j+1)*nzz + k*nxx*nzz]) / dx +
                               (bp1*adjoint_comp[i + j*nzz + (k-1)*nxx*nzz] - bm2*adjoint_comp[i + j*nzz + (k+1)*nxx*nzz]) / dy +
                               (cp1*adjoint_comp[(i-1) + j*nzz + k*nxx*nzz] - cm2*adjoint_comp[(i+1) + j*nzz + k*nxx*nzz]) / dz;

                    float fg = (eg + source_grad[i + j*nzz + k*nxx*nzz]) / d;
                    float fc = (ec + source_comp[i + j*nzz + k*nxx*nzz]) / d;

                    if (adjoint_grad[i + j*nzz + k*nxx*nzz] > fg) 
                        adjoint_grad[i + j*nzz + k*nxx*nzz] = fg;

                    if (adjoint_comp[i + j*nzz + k*nxx*nzz] > fc) 
                        adjoint_comp[i + j*nzz + k*nxx*nzz] = fc;
                }
            }
        }
    }
}
