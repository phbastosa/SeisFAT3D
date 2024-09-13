# include "kirchhoff.cuh"
# include <omp.h>

void Kirchhoff::set_modeling_type()
{
    modeling = new Ultimate_FSM();
}

void Kirchhoff::image_building()
{
    modeling->set_runtime();

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

        for (int node = 0; node < modeling->total_nodes; node++)
        {            
            import_binary_float("../outputs/snapshots/travel_time_volume_receiver_" + std::to_string(node+1) + ".bin", Tr, modeling->nPoints);

            #pragma omp parallel for     
            for (int index = 0; index < modeling->nPoints; index++)
            {   
                int k = (int) (index / (modeling->nx*modeling->nz));                  // y direction
                int j = (int) (index - k*modeling->nx*modeling->nz) / modeling->nz;   // x direction
                int i = (int) (index - j*modeling->nz - k*modeling->nx*modeling->nz); // z direction                

                Im[index] = modeling->wavefield_output[index] + Tr[index] - modeling->receiver_output[node];

                int seism_index = (int)(Im[index] / dt);

                if (seism_index < nt) image[index] += data[seism_index + node*nt];
            }
        }
    }
 
    export_binary_float("../outputs/images/image_kirchhoff_migration_test.bin", image, modeling->nPoints);

    modeling->get_runtime();
}