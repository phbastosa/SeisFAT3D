# include "ADKDM.cuh"

void ADKDM::set_migration()
{
    domain = "Angle Domain";
    migType = "ADKDM";
    m_samples = modeling->nz*nang*nCMP;

    output_path = images_folder + migType + "_result_" + std::to_string(modeling->nz) + "x" + std::to_string(nCMP) + "x" + std::to_string(nang) + ".bin";
}

void ADKDM::perform_forward()
{
    angle_domain_forward_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_model, modeling->dx, modeling->dz, dt, da, modeling->nxx, modeling->nzz, nt, nang, modeling->nb, cmpId);
}

void ADKDM::perform_adjoint()
{
    angle_domain_adjoint_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_model, modeling->dx, modeling->dz, dt, da, modeling->nxx, modeling->nzz, nt, nang, modeling->nb, cmpId);
}
