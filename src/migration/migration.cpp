# include "migration.hpp"

void Migration::set_parameters()
{
    set_modeling_type();

    modeling->file = file;

    modeling->set_parameters();





    modeling->set_runtime();
} 

void Migration::import_gathers()
{



}

void Migration::print_information()
{
    auto clear = system("clear");

    std::cout << "\033[1mSeisFAT3D\033[m - Migration program\n\n";

    std::cout << "Model dimensions: (z = " << (modeling->nz-1)*modeling->dz << 
                                  ", x = " << (modeling->nx-1)*modeling->dx << 
                                  ", y = " << (modeling->ny-1)*modeling->dy << ") m\n\n";

    std::cout << "Migration type: \033[1m" << type_message << "\033[m\n\n";

    std::cout << "Running shot " << modeling->shot_index+1 << " of " << modeling->total_shots;

    std::cout << " at position (z = " << modeling->geometry->shots.z[modeling->shot_index] << 
                             ", x = " << modeling->geometry->shots.x[modeling->shot_index] << 
                             ", y = " << modeling->geometry->shots.y[modeling->shot_index] << ") m\n\n";
}

void Migration::cross_correlation()
{
    for (int shot = 0; shot < modeling->total_shots; shot++)
    {
        modeling->shot_index = shot;

        print_information();

        source_to_receiver_propagation();
        receiver_to_source_propagation();
    }
}

void Migration::source_to_receiver_propagation()
{
    modeling->set_initial_conditions();

    modeling->forward_propagation();
        

}

void Migration::receiver_to_source_propagation()
{


}

void Migration::image_compensation()
{




}

void Migration::export_seismic_volume()
{
    output_image_file = output_image_folder + "image_etc";

    export_binary_float(output_image_file, image, modeling->nPoints);

    modeling->get_runtime();
}