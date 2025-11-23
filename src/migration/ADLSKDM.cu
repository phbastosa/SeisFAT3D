# include "ADLSKDM.cuh"

void ADLSKDM::set_migration()
{
    domain = "Angle Domain";
    migType = "ADLSKDM";
    m_samples = old_nz*nCMP*nang;

    output_path = images_folder + migType + "_result_" + std::to_string(old_nz) + "x" + std::to_string(nCMPy) + "x" + std::to_string(nCMPx) + "x" + std::to_string(nang) + "_iteration_" + std::to_string(max_it) + ".bin";
}

void ADLSKDM::perform_forward()
{
    angle_domain_forward_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_model, modeling->dx, modeling->dz, dt, da, modeling->nxx, modeling->nzz, nt, nang, modeling->nb, cmpId);
}

void ADLSKDM::perform_adjoint()
{
    angle_domain_adjoint_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_model, modeling->dx, modeling->dz, dt, da, modeling->nxx, modeling->nzz, nt, nang, modeling->nb, cmpId);
}

void ADLSKDM::perform_adjoint_gradient()
{
    angle_domain_adjoint_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_gradient, modeling->dx, modeling->dz, dt, da, modeling->nxx, modeling->nzz, nt, nang, modeling->nb, cmpId);
}

void ADLSKDM::perform_forward_direction()
{
    angle_domain_forward_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_direction, modeling->dx, modeling->dz, dt, da, modeling->nxx, modeling->nzz, nt, nang, modeling->nb, cmpId);    
}
