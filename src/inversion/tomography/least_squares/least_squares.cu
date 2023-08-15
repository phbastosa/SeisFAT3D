# include "least_squares.cuh"

void Least_Squares::set_parameters()
{
    set_general_parameters();
    
    set_forward_modeling();

    set_main_components();

    dx_tomo = std::stof(catch_parameter("dx_tomo", file));
    dy_tomo = std::stof(catch_parameter("dy_tomo", file));
    dz_tomo = std::stof(catch_parameter("dz_tomo", file));

    nz_tomo = (int)((modeling->nz-1) * modeling->dz / dz_tomo) + 1;    
    nx_tomo = (int)((modeling->nx-1) * modeling->dx / dx_tomo) + 1;    
    ny_tomo = (int)((modeling->ny-1) * modeling->dy / dy_tomo) + 1;  

    n_model = nx_tomo * ny_tomo * nz_tomo;

    std::cout<<"Least squares first arrival tomography"<<std::endl;
}

void Least_Squares::forward_modeling()
{
    for (int shot = 0; shot < modeling->total_shots; shot++)
    {
        modeling->shot_id = shot;

        modeling->info_message();

        tomography_message();

        modeling->initial_setup();
        modeling->forward_solver();
        modeling->build_outputs();

        extract_calculated_data();

        // gradient_ray_tracing()
    }
}

void Least_Squares::gradient_ray_tracing()
{




}

void Least_Squares::optimization()
{


}

