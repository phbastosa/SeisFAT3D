# include "IDKDM.cuh"

void IDKDM::set_migration()
{
    domain = "Image Domain";
    migType = "IDKDM";
    m_samples = old_nPoints;

    output_path = seismic_folder + migType + "_result_" + std::to_string(old_nz) + "x" + std::to_string(old_nx) + "x" + std::to_string(old_ny) + ".bin";
}

void IDKDM::perform_adjoint()
{
    const int TILE = 8;
    const int HALO = 4;

    dim3 blk(TILE, TILE, TILE);

    dim3 grd((old_nz + blk.x - 1) / blk.x,
             (old_nx + blk.y - 1) / blk.y,
             (old_ny + blk.z - 1) / blk.z);

    const float tile_phys_z = blk.x * old_dh;
    const float tile_phys_x = blk.y * old_dh;
    const float tile_phys_y = blk.z * old_dh;

    const int sm_z = (int)ceilf(tile_phys_z / new_dh) + HALO + 1;
    const int sm_x = (int)ceilf(tile_phys_x / new_dh) + HALO + 1;
    const int sm_y = (int)ceilf(tile_phys_y / new_dh) + HALO + 1;

    const int total_sm_nodes = sm_z * sm_x * sm_y;

    const int sm_bytes = 3 * total_sm_nodes * sizeof(float);

    image_domain_adjoint_kernel<<<grd,blk,sm_bytes,0>>>(
        modeling->d_S, d_Ts, d_Tr, d_data, d_model, d_xsrc, d_ysrc, 
        d_xrec, d_yrec, aperture, max_offset, sm_z, sm_x, sm_y, new_dh, 
        old_dh, dt, modeling->nxx, modeling->nyy, modeling->nzz, 
        modeling->nb, old_nx, old_ny, old_nz, nt, nsy, nrx);
}

void IDKDM::perform_forward()
{
    const int TILE = 8;
    const int HALO = 4;

    dim3 blk(TILE, TILE, TILE);

    dim3 grd((old_nz + blk.x - 1) / blk.x,
             (old_nx + blk.y - 1) / blk.y,
             (old_ny + blk.z - 1) / blk.z);

    const float tile_phys_z = blk.x * old_dh;
    const float tile_phys_x = blk.y * old_dh;
    const float tile_phys_y = blk.z * old_dh;

    const int sm_z = (int)ceilf(tile_phys_z / new_dh) + HALO + 1;
    const int sm_x = (int)ceilf(tile_phys_x / new_dh) + HALO + 1;
    const int sm_y = (int)ceilf(tile_phys_y / new_dh) + HALO + 1;

    const int total_sm_nodes = sm_z * sm_x * sm_y;

    const int sm_bytes = 3 * total_sm_nodes * sizeof(float);

    image_domain_forward_kernel<<<grd,blk,sm_bytes,0>>>(
        modeling->d_S, d_Ts, d_Tr, d_data, d_model, d_xsrc, d_ysrc, 
        d_xrec, d_yrec, aperture, max_offset, sm_z, sm_x, sm_y, new_dh, 
        old_dh, dt, modeling->nxx, modeling->nyy, modeling->nzz, 
        modeling->nb, old_nx, old_ny, old_nz, nt, nsy, nrx);
}
