# include "IDLSKDM.cuh"

void IDLSKDM::set_migration()
{
    domain = "Image Domain";
    migType = "IDLSKDM";
    m_samples = old_nPoints;

    output_path = seismic_folder + migType + "_result_" + std::to_string(old_nz) + "x" + std::to_string(old_nx) + "x" + std::to_string(old_ny) + "_iteration_" + std::to_string(max_it) + ".bin";
}

void IDLSKDM::perform_adjoint()
{
    image_domain_adjoint_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_model, dt, nt, old_dx, old_dy, old_dz, new_dx, new_dy, 
                                                      new_dz, old_nx, old_ny, old_nz, modeling->nxx, modeling->nyy, modeling->nzz, modeling->nb,
                                                      aperture, CMPx, CMPy);
}

void IDLSKDM::perform_forward()
{
    image_domain_forward_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_model, dt, nt, old_dx, old_dy, old_dz, new_dx, new_dy, 
                                                      new_dz, old_nx, old_ny, old_nz, modeling->nxx, modeling->nyy, modeling->nzz, modeling->nb,
                                                      aperture, CMPx, CMPy);
}

void IDLSKDM::perform_adjoint_gradient()
{
    image_domain_adjoint_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_gradient, dt, nt, old_dx, old_dy, old_dz, new_dx, new_dy, 
                                                      new_dz, old_nx, old_ny, old_nz, modeling->nxx, modeling->nyy, modeling->nzz, modeling->nb,
                                                      aperture, CMPx, CMPy);
}

void IDLSKDM::perform_forward_direction()
{
    image_domain_forward_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_direction, dt, nt, old_dx, old_dy, old_dz, new_dx, new_dy, 
                                                      new_dz, old_nx, old_ny, old_nz, modeling->nxx, modeling->nyy, modeling->nzz, modeling->nb,
                                                      aperture, CMPx, CMPy);    
}
