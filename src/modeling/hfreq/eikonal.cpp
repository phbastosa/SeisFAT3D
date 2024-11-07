# include "eikonal.hpp"

void Eikonal::set_specifications()
{
    set_properties();    
    set_conditions();    

    synthetic_data = new float[max_spread]();
}

void Eikonal::set_boundaries()
{
    nb = 1;    
}

void Eikonal::initialization()
{
    int sId = geometry->sInd[srcId];

    sIdx = (int)(geometry->xsrc[sId] / dx) + nb;
    sIdy = (int)(geometry->ysrc[sId] / dy) + nb;
    sIdz = (int)(geometry->zsrc[sId] / dz) + nb;

    # pragma omp parallel for
    for (int index = 0; index < volsize; index++) 
        T[index] = 1e6f;

    for (int k = 0; k < 3; k++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int i = 0; i < 3; i++)
            {
                int yi = sIdy + (k - 1);
                int xi = sIdx + (j - 1);
                int zi = sIdz + (i - 1);

                T[zi + xi*nzz + yi*nxx*nzz] = S[zi + xi*nzz + yi*nxx*nzz] * 
                    sqrtf(powf((xi - nb)*dx - geometry->xsrc[sId], 2.0f) + 
                          powf((yi - nb)*dz - geometry->ysrc[sId], 2.0f) +
                          powf((zi - nb)*dz - geometry->zsrc[sId], 2.0f));
            }
        }
    }
}

void Eikonal::compute_seismogram()
{
    int spread = 0;

    int sId = geometry->sInd[srcId];

    for (recId = geometry->iRec[sId]; recId < geometry->fRec[sId]; recId++)
    {
        float x = geometry->xrec[recId];
        float y = geometry->yrec[recId];
        float z = geometry->zrec[recId];

        float x0 = floorf(x / dx) * dx;
        float y0 = floorf(y / dy) * dy;
        float z0 = floorf(z / dz) * dz;

        float x1 = floorf(x / dx) * dx + dx;
        float y1 = floorf(y / dy) * dy + dy;
        float z1 = floorf(z / dz) * dz + dz;

        int id = ((int)(z / dz) + nb) + ((int)(x / dx) + nb)*nzz + ((int)(y / dy) + nb)*nxx*nzz;

        float c000 = T[id];
        float c001 = T[id + 1];
        float c100 = T[id + nzz]; 
        float c101 = T[id + 1 + nzz]; 
        float c010 = T[id + nxx*nzz]; 
        float c011 = T[id + 1 + nxx*nzz]; 
        float c110 = T[id + nzz + nxx*nzz]; 
        float c111 = T[id + 1 + nzz + nxx*nzz];

        float xd = (x - x0) / (x1 - x0);
        float yd = (y - y0) / (y1 - y0);
        float zd = (z - z0) / (z1 - z0);

        float c00 = c000*(1 - xd) + c100*xd;    
        float c01 = c001*(1 - xd) + c101*xd;    
        float c10 = c010*(1 - xd) + c110*xd;    
        float c11 = c011*(1 - xd) + c111*xd;    

        float c0 = c00*(1 - yd) + c10*yd;
        float c1 = c01*(1 - yd) + c11*yd;

        synthetic_data[spread++] = c0*(1 - zd) + c1*zd;
    }
}

void Eikonal::export_synthetic_data()
{
    int sId = geometry->sInd[srcId];    
    std::string data_file = data_folder + modeling_type + "nStations" + std::to_string(geometry->spread[sId]) + "_shot_" + std::to_string(sId+1) + ".bin";
    export_binary_float(data_file, synthetic_data, geometry->spread[sId]);    
}