# include "eikonal.cuh"

void Eikonal::FSM_parameters()
{
    padb = 1;

    title = "Eikonal solver for acoustic isotropic media\n\nSolving eikonal equation with the \033[32mNoble, Gesret and Belayouni (2014)\033[0;0m formulation\n";    
}

void Eikonal::FSM_components()
{
    nxx = nx + 2*padb;
    nyy = ny + 2*padb;
    nzz = nz + 2*padb;
}

void Eikonal::FSM_init()
{
    int sId = sidz + sidx*nzz + sidy*nxx*nzz;

    for (int index = 0; index < volsize; index++)
        T[index] = 1e6f;

    // Neighboring source points initialization with analitical traveltime

    T[sId] = S[sId] * sqrtf(powf((sidx-padb)*dh - geometry->shots.x[shot_id], 2.0f) + powf((sidy-padb)*dh - geometry->shots.y[shot_id], 2.0f) + powf((sidz-padb)*dh - geometry->shots.z[shot_id], 2.0f));
    
    T[sId + 1] = S[sId] * sqrtf(powf((sidx-padb)*dh - geometry->shots.x[shot_id], 2.0f) + powf((sidy-padb)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)+1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - 1] = S[sId] * sqrtf(powf((sidx-padb)*dh - geometry->shots.x[shot_id], 2.0f) + powf((sidy-padb)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)-1)*dh - geometry->shots.z[shot_id], 2.0f));

    T[sId + nzz] = S[sId] * sqrtf(powf(((sidx-padb)+1)*dh - geometry->shots.x[shot_id], 2.0f) + powf((sidy-padb)*dh - geometry->shots.y[shot_id], 2.0f) + powf((sidz-padb)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - nzz] = S[sId] * sqrtf(powf(((sidx-padb)-1)*dh - geometry->shots.x[shot_id], 2.0f) + powf((sidy-padb)*dh - geometry->shots.y[shot_id], 2.0f) + powf((sidz-padb)*dh - geometry->shots.z[shot_id], 2.0f));
    
    T[sId + nxx*nzz] = S[sId] * sqrtf(powf((sidx-padb)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)+1)*dh - geometry->shots.y[shot_id], 2.0f) + powf((sidz-padb)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - nxx*nzz] = S[sId] * sqrtf(powf((sidx-padb)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)-1)*dh - geometry->shots.y[shot_id], 2.0f) + powf((sidz-padb)*dh - geometry->shots.z[shot_id], 2.0f));
    
    T[sId + 1 + nzz] = S[sId] * sqrtf(powf(((sidx-padb)+1)*dh - geometry->shots.x[shot_id], 2.0f) + powf((sidy-padb)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)+1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId + 1 - nzz] = S[sId] * sqrtf(powf(((sidx-padb)+1)*dh - geometry->shots.x[shot_id], 2.0f) + powf((sidy-padb)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)-1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - 1 + nzz] = S[sId] * sqrtf(powf(((sidx-padb)-1)*dh - geometry->shots.x[shot_id], 2.0f) + powf((sidy-padb)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)+1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - 1 - nzz] = S[sId] * sqrtf(powf(((sidx-padb)-1)*dh - geometry->shots.x[shot_id], 2.0f) + powf((sidy-padb)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)-1)*dh - geometry->shots.z[shot_id], 2.0f));
    
    T[sId + 1 + nxx*nzz] = S[sId] * sqrtf(powf((sidx-padb)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)+1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)+1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId + 1 - nxx*nzz] = S[sId] * sqrtf(powf((sidx-padb)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)-1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)+1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - 1 + nxx*nzz] = S[sId] * sqrtf(powf((sidx-padb)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)+1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)-1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - 1 - nxx*nzz] = S[sId] * sqrtf(powf((sidx-padb)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)-1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)-1)*dh - geometry->shots.z[shot_id], 2.0f));
    
    T[sId + nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)+1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)+1)*dh - geometry->shots.y[shot_id], 2.0f) + powf((sidz-padb)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId + nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)+1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)-1)*dh - geometry->shots.y[shot_id], 2.0f) + powf((sidz-padb)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)-1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)+1)*dh - geometry->shots.y[shot_id], 2.0f) + powf((sidz-padb)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)-1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)-1)*dh - geometry->shots.y[shot_id], 2.0f) + powf((sidz-padb)*dh - geometry->shots.z[shot_id], 2.0f));
    
    T[sId + 1 + nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)+1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)+1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)+1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId + 1 - nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)-1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)+1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)+1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId + 1 + nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)+1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)-1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)+1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId + 1 - nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)-1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)-1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)+1)*dh - geometry->shots.z[shot_id], 2.0f));

    T[sId - 1 + nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)+1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)+1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)-1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - 1 - nzz + nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)-1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)+1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)-1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - 1 + nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)+1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)-1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)-1)*dh - geometry->shots.z[shot_id], 2.0f));
    T[sId - 1 - nzz - nxx*nzz] = S[sId] * sqrtf(powf(((sidx-padb)-1)*dh - geometry->shots.x[shot_id], 2.0f) + powf(((sidy-padb)-1)*dh - geometry->shots.y[shot_id], 2.0f) + powf(((sidz-padb)-1)*dh - geometry->shots.z[shot_id], 2.0f));

    dzi = 1.0f / dh;
    dxi = 1.0f / dh;
    dyi = 1.0f / dh;

    dz2i = 1.0f / (dh*dh);
    dx2i = 1.0f / (dh*dh);
    dy2i = 1.0f / (dh*dh);

    dz2dx2 = dz2i * dx2i;
    dz2dy2 = dz2i * dy2i;
    dx2dy2 = dx2i * dy2i;

    dsum = dz2i + dx2i + dy2i;
}

void Eikonal::FSM_solver()
{
    init_sweep();
    full_sweep();
}

void Eikonal::init_sweep()
{
    // First sweeping: Top->Bottom; West->East; South->North
    sgntz = 1; sgntx = 1; sgnty = 1; 
    sgnvz = 1; sgnvx = 1; sgnvy = 1;

    for (k = max(1, sidy); k < nyy; k++)
    {
        for (j = max(1, sidx); j < nxx; j++)
        {
            for (i = max(1, sidz); i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Second sweeping: Top->Bottom; East->West; South->North
    sgntz = -1; sgntx = 1; sgnty = 1;
    sgnvz =  0; sgnvx = 1; sgnvy = 1;

    for (k = max(1, sidy); k < nyy; k++)
    {
        for (j = max(1, sidx); j < nxx; j++)
        {
            for (i = sidz + 1; i >= 0 ; i--)
            {
                inner_sweep();
            }
        }
    }
    
    // Third sweeping: Top->Bottom; West->East; North->South
    sgntz = 1; sgntx = 1; sgnty = -1;
    sgnvz = 1; sgnvx = 1; sgnvy =  0;

    for (k = sidy + 1; k >= 0; k--)
    {
        for (j = max(1, sidx); j < nxx; j++)
        {
            for (i = max(1, sidz); i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Fourth sweeping: Top->Bottom ; East->West ; North->South
    sgntz = -1; sgntx = 1; sgnty = -1;
    sgnvz =  0; sgnvx = 1; sgnvy =  0;

    for (k = sidy + 1; k >= 0; k--)
    {
        for (j = max(1, sidx); j < nxx; j++)
        {
            for (i = sidz + 1; i >= 0 ; i--)
            {
                inner_sweep();
            }
        }
    }

    // Fifth sweeping: Bottom->Top; West->East; South->North
    sgntz = 1; sgntx = -1; sgnty = 1;
    sgnvz = 1; sgnvx =  0; sgnvy = 1;

    for (k = max(1, sidy); k < nyy; k++)
    {
        for (j = sidx + 1; j >= 0; j--)
        {
            for (i = max(1, sidz); i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Sixth sweeping: Bottom->Top; East->West; South->North
    sgntz = -1; sgntx = -1; sgnty = 1;
    sgnvz =  0; sgnvx =  0; sgnvy = 1;

    for (k = max(1, sidy); k < nyy; k++)
    {
        for (j = sidx + 1; j >= 0; j--)
        {
            for (i = sidz + 1; i >= 0; i--)
            {
                inner_sweep();
            }
        }
    }

    // Seventh sweeping: Bottom->Top; West->East; North->South
    sgntz = 1; sgntx = -1; sgnty = -1;
    sgnvz = 1; sgnvx =  0; sgnvy =  0;

    for (k = sidy + 1; k >= 0; k--)
    {
        for (j = sidx + 1; j >= 0; j--)
        {
            for (i = max(1, sidz); i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Eighth sweeping: Bottom->Top; East->West; North->South
    sgntz = -1; sgntx = -1; sgnty = -1;
    sgnvz =  0; sgnvx =  0; sgnvy =  0;

    for (k = sidy + 1; k >= 0; k--)
    {
        for (j = sidx + 1; j >= 0; j--)
        {
            for (i = sidz + 1; i >= 0; i--)
            {
                inner_sweep();
            }
        }
    }
}

void Eikonal::full_sweep()
{
    // First sweeping: Top->Bottom; West->East; South->North 
    sgntz = 1; sgntx = 1; sgnty = 1; 
    sgnvz = 1; sgnvx = 1; sgnvy = 1;

    for (k = 1; k < nyy; k++)
    {
        for (j = 1; j < nxx; j++)
        {
            for (i = 1; i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Second sweeping: Top->Bottom; East->West; South->North
    sgntz = -1; sgntx = 1; sgnty = 1;
    sgnvz =  0; sgnvx = 1; sgnvy = 1;

    for (k = 1; k < nyy; k++)
    {
        for (j = 1; j < nxx; j++)
        {
            for (i = nzz - 2; i >= 0; i--)
            {
                inner_sweep();
            }
        }
    }
    
    // Third sweeping: Top->Bottom; West->East; North->South
    sgntz = 1; sgntx = 1; sgnty = -1;
    sgnvz = 1; sgnvx = 1; sgnvy =  0;

    for (k = nyy - 2; k >= 0; k--)
    {
        for (j = 1; j < nxx; j++)
        {
            for (i = 1; i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Fourth sweeping: Top->Bottom ; East->West ; North->South
    sgntz = -1; sgntx = 1; sgnty = -1;
    sgnvz =  0; sgnvx = 1; sgnvy =  0;

    for (k = nyy - 2; k >= 0; k--)
    {
        for (j = 1; j < nxx; j++)
        {
            for (i = nzz - 2; i >= 0; i--)
            {
                inner_sweep();
            }
        }
    }

    // Fifth sweeping: Bottom->Top; West->East; South->North
    sgntz = 1; sgntx = -1; sgnty = 1;
    sgnvz = 1; sgnvx =  0; sgnvy = 1;

    for (k = 1; k < nyy; k++)
    {
        for (j = nxx - 2; j >= 0; j--)
        {
            for (i = 1; i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Sixth sweeping: Bottom->Top; East->West; South->North
    sgntz = -1; sgntx = -1; sgnty = 1;
    sgnvz =  0; sgnvx =  0; sgnvy = 1;

    for (k = 1; k < nyy; k++)
    {
        for (j = nxx - 2; j >= 0; j--)
        {
            for (i = nzz - 2; i >= 0; i--)
            {
                inner_sweep();
            }
        }
    }

    // Seventh sweeping: Bottom->Top; West->East; North->South
    sgntz = 1; sgntx = -1; sgnty = -1;
    sgnvz = 1; sgnvx =  0; sgnvy =  0;

    for (k = nyy - 2; k >= 0; k--)
    {
        for (j = nxx - 2; j >= 0; j--)
        {
            for (i = 1; i < nzz; i++)
            {
                inner_sweep();
            }
        }
    }

    // Eighth sweeping: Bottom->Top; East->West; North->South
    sgntz = -1; sgntx = -1; sgnty = -1;
    sgnvz =  0; sgnvx =  0; sgnvy =  0;

    for (k = nyy - 2; k >= 0; k--)
    {
        for (j = nxx - 2; j >= 0; j--)
        {
            for (i = nzz - 2; i >= 0; i--)
            {
                inner_sweep();
            }
        }
    }
}

void Eikonal::inner_sweep()
{
    float ta, tb, tc, t1, t2, t3, Sref;
    float t1D1, t1D2, t1D3, t1D, t2D1, t2D2, t2D3, t2D, t3D;

    // Index of velocity nodes
    int i1 = i - sgnvz; 
    int j1 = j - sgnvx; 
    int k1 = k - sgnvy;

    // Get local times of surrounding points
    float tv = T[(i - sgntz) + j*nzz + k*nxx*nzz];
    float te = T[i + (j - sgntx)*nzz + k*nxx*nzz];
    float tn = T[i + j*nzz + (k - sgnty)*nxx*nzz];
    float tev = T[(i - sgntz) + (j - sgntx)*nzz + k*nxx*nzz];
    float ten = T[i + (j - sgntx)*nzz + (k - sgnty)*nxx*nzz];
    float tnv = T[(i - sgntz) + j*nzz + (k - sgnty)*nxx*nzz];
    float tnve = T[(i - sgntz) + (j - sgntx)*nzz + (k - sgnty)*nxx*nzz];     

    int ijk = i + j*nzz + k*nxx*nzz;

    //------------------- 1D operators ---------------------------------------------------------------------------------------------------
    t1D1 = 1e5; t1D2 = 1e5; t1D3 = 1e5;     

    // Z direction
    t1D1 = tv + dh * min(S[i1 + max(j-1,1)*nzz   + max(k-1,1)*nxx*nzz], 
                     min(S[i1 + max(j-1,1)*nzz   + min(k,nyy-1)*nxx*nzz], 
                     min(S[i1 + min(j,nxx-1)*nzz + max(k-1,1)*nxx*nzz],
                         S[i1 + min(j,nxx-1)*nzz + min(k,nyy-1)*nxx*nzz]))); 

    // X direction
    t1D2 = te + dh * min(S[max(i-1,1)   + j1*nzz + max(k-1,1)*nxx*nzz], 
                     min(S[min(i,nzz-1) + j1*nzz + max(k-1,1)*nxx*nzz],
                     min(S[max(i-1,1)   + j1*nzz + min(k,nyy-1)*nxx*nzz], 
                         S[min(i,nzz-1) + j1*nzz + min(k,nyy-1)*nxx*nzz])));

    // Y direction
    t1D3 = tn + dh * min(S[max(i-1,1)   + max(j-1,1)*nzz   + k1*nxx*nzz], 
                     min(S[max(i-1,1)   + min(j,nxx-1)*nzz + k1*nxx*nzz],
                     min(S[min(i,nzz-1) + max(j-1,1)*nzz   + k1*nxx*nzz], 
                         S[min(i,nzz-1) + min(j,nxx-1)*nzz + k1*nxx*nzz])));

    t1D = min(t1D1, min(t1D2, t1D3));

    //------------------- 2D operators - 4 points operator ---------------------------------------------------------------------------------------------------
    t2D1 = 1e6f; t2D2 = 1e6f; t2D3 = 1e6f;

    // XZ plane ----------------------------------------------------------------------------------------------------------------------------------------------
    Sref = min(S[i1 + j1*nzz + max(k-1,1)*nxx*nzz], S[i1 + j1*nzz + min(k, nyy-1)*nxx*nzz]);
    
    if ((tv < te + dh*Sref) && (te < tv + dh*Sref))
    {
        ta = tev + te - tv;
        tb = tev - te + tv;

        t2D1 = ((tb*dz2i + ta*dx2i) + sqrtf(4.0f*Sref*Sref*(dz2i + dx2i) - dz2i*dx2i*(ta - tb)*(ta - tb))) / (dz2i + dx2i);
    }

    // YZ plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
    Sref = min(S[i1 + max(j-1,1)*nzz + k1*nxx*nzz], S[i1 + min(j,nxx-1)*nzz + k1*nxx*nzz]);

    if((tv < tn + dh*Sref) && (tn < tv + dh*Sref))
    {
        ta = tv - tn + tnv;
        tb = tn - tv + tnv;
        
        t2D2 = ((ta*dz2i + tb*dy2i) + sqrtf(4.0f*Sref*Sref*(dz2i + dy2i) - dz2i*dy2i*(ta - tb)*(ta - tb))) / (dz2i + dy2i); 
    }

    // XY plane -------------------------------------------------------------------------------------------------------------------------------------------------------------
    Sref = min(S[max(i-1,1) + j1*nzz + k1*nxx*nzz], S[min(i,nzz-1) + j1*nzz + k1*nxx*nzz]);

    if((te < tn + dh*Sref) && (tn < te + dh*Sref))
    {
        ta = te - tn + ten;
        tb = tn - te + ten;

        t2D3 = ((ta*dx2i + tb*dy2i) + sqrtf(4.0f*Sref*Sref*(dx2i + dy2i) - dx2i*dy2i*(ta - tb)*(ta - tb))) / (dx2i + dy2i);
    }

    t2D = min(t2D1, min(t2D2, t2D3));

    //------------------- 3D operators ---------------------------------------------------------------------------------------------------
    t3D = 1e6f;

    Sref = S[i1 + j1*nzz + k1*nxx*nzz];

    ta = te - 0.5f*tn + 0.5f*ten - 0.5f*tv + 0.5f*tev - tnv + tnve;
    tb = tv - 0.5f*tn + 0.5f*tnv - 0.5f*te + 0.5f*tev - ten + tnve;
    tc = tn - 0.5f*te + 0.5f*ten - 0.5f*tv + 0.5f*tnv - tev + tnve;

    if (min(t1D,t2D) > max(tv, max(te, tn)))
    {
        t2 = 9.0f*Sref*Sref*dsum;
        
        t3 = dz2dx2*(ta - tb)*(ta - tb) + dz2dy2*(tb - tc)*(tb - tc) + dx2dy2*(ta - tc)*(ta - tc);
        
        if (t2 >= t3)
        {
            t1 = tb*dz2i + ta*dx2i + tc*dy2i;        
            
            t3D = (t1 + sqrtf(t2 - t3)) / dsum;
        }
    }
   
    T[ijk] = min(T[ijk], min(t1D, min(t2D, t3D)));
}
