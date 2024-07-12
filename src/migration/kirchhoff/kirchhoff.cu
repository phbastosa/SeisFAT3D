# include "kirchhoff.cuh"

void Kirchhoff::set_modeling_type()
{
    modeling = new Ultimate_FSM();
}

void Kirchhoff::image_building()
{
    modeling->set_runtime();

    float pi = 4.0f*atanf(1.0f);

    // for (int node = 0; node < modeling->total_nodes; node++)
    // {            
    //     modeling->get_information();

    //     std::cout << "Computing travel times at receiver " << node+1 << " of " << modeling->total_nodes << "\n";

    //     modeling->sidx = (int)(modeling->geometry->nodes.x[node] / modeling->dx) + modeling->nbxl;
    //     modeling->sidy = (int)(modeling->geometry->nodes.y[node] / modeling->dy) + modeling->nbyl;
    //     modeling->sidz = (int)(modeling->geometry->nodes.z[node] / modeling->dz) + modeling->nbzu; 

    //     modeling->source_index = modeling->sidz + modeling->sidx*modeling->nzz + modeling->sidy*modeling->nxx*modeling->nzz;

    //     modeling->initialization();
    //     modeling->set_forward_solver();

    //     export_binary_float("../outputs/snapshots/travel_time_volume_receiver_" + std::to_string(node+1) + ".bin", modeling->wavefield_output, modeling->nPoints);
    // }

    for (int shot = 0; shot < 1; shot++)
    {
        modeling->shot_index = shot;

        read_input_data();
        
        modeling->get_information();

        std::cout << "Kirchhoff depth migration \n\n";

        modeling->set_configuration();
        modeling->set_forward_solver();

        for (int node = 0; node < 1; node++)
        {            
            import_binary_float("../outputs/snapshots/travel_time_volume_receiver_" + std::to_string(node+1) + ".bin", Tr, modeling->nPoints);
        
            for (int index = 0; index < modeling->nPoints; index++)
            {   
                int k = (int) (index / (modeling->nx*modeling->nz));                  // y direction
                int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;   // x direction
                int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz); // z direction                

                Im[index] = modeling->wavefield_output[index] + Tr[index] - modeling->receiver_output[node];

                int seism_index = (int)(Im[index] / dt);

                if ((i > 0) && (i < modeling->nz-1) && (j > 0) && (j < modeling->nx-1) && (k > 0) && (k < modeling->ny-1))
                {
                    float dTs_dz = (modeling->wavefield_output[index + 1] - Ts[index - 1]) / (2.0f*modeling->dz);     
                    float dTs_dx = (modeling->wavefield_output[index + modeling->nz] - Ts[index - modeling->nz]) / (2.0f*modeling->dx);
                    float dTs_dy = (modeling->wavefield_output[index + modeling->nx*modeling->nz] - Ts[index - modeling->nx*modeling->nz]) / (2.0f*modeling->dy);

                    float dTs_norm = sqrtf(dTs_dx*dTs_dx + dTs_dy*dTs_dy + dTs_dz*dTs_dz);

                    float dTr_dz = (Tr[index + 1] - Tr[index - 1]) / (2.0f*modeling->dz);
                    float dTr_dx = (Tr[index + modeling->nz] - Tr[index - modeling->nz]) / (2.0f*modeling->dx);
                    float dTr_dy = (Tr[index + modeling->nx*modeling->nz] - Tr[index - modeling->nx*modeling->nz]) / (2.0f*modeling->dy);

                    float dTr_norm = sqrtf(dTr_dx*dTr_dx + dTr_dy*dTr_dy + dTr_dz*dTr_dz);

                    dTx[index] = dTs_dx/dTs_norm + dTr_dx/dTr_norm; 
                    dTy[index] = dTs_dy/dTs_norm + dTr_dy/dTr_norm; 
                    dTz[index] = dTs_dz/dTs_norm + dTr_dz/dTr_norm; 
                        
                    //     float Ts_dot_Tr = dTs_dx*dTr_dx + dTs_dy*dTr_dy + dTs_dz*dTr_dz;
                    //     float norm_dots = sqrtf(dTs_dx*dTs_dx + dTs_dy*dTs_dy + dTs_dz*dTs_dz) * sqrtf(dTr_dx*dTr_dx + dTr_dy*dTr_dy + dTr_dz*dTr_dz); 

                    //     float angle = 0.5f*acosf(Ts_dot_Tr / norm_dots);
                }
            
                if ((seism_index < nt)) 
                    image[index] += data[seism_index + node*nt];  
            }

        }

        export_binary_float("test_image.bin", image, modeling->nPoints);
        
        export_binary_float("test_dTx.bin", dTx, modeling->nPoints);
        export_binary_float("test_dTy.bin", dTy, modeling->nPoints);
        export_binary_float("test_dTz.bin", dTz, modeling->nPoints);
    }

    modeling->get_runtime();
}