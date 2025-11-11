# include "IDKDM.cuh"

void IDKDM::set_migration()
{
    domain = "Image Domain";
    migType = "IDKDM";
    m_samples = modeling->nPoints;

    output_path = images_folder + migType + "_result_" + std::to_string(modeling->nz) + "x" + std::to_string(modeling->nx) + ".bin";
}

void IDKDM::perform_forward()
{
    image_domain_forward_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_model, modeling->dx, modeling->dz, dt, modeling->nxx, modeling->nzz, nt, modeling->nb);
}

void IDKDM::perform_adjoint()
{
    image_domain_adjoint_kernel<<<nBlocks,NTHREADS>>>(modeling->d_S, d_Ts, d_Tr, d_data, d_model, modeling->dx, modeling->dz, dt, modeling->nxx, modeling->nzz, nt, modeling->nb);
}
