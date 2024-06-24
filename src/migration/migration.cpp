# include "migration.hpp"

void Migration::set_parameters()
{
    set_modeling_type();

    modeling->file = file;

    modeling->set_parameters();

    total_shots = modeling->total_shots;
    total_nodes = modeling->total_nodes;

    modeling->set_runtime();





} 

void Migration::import_seismic_data()
{



}

void Migration::print_information()
{
    modeling->print_information();

    std::cout<<"testing new functionalities\n";



}

void Migration::source_to_receiver_propagation()
{
    modeling->shot_index = shot_index;

    print_information();
    
    modeling->set_initial_conditions();

    modeling->forward_propagation();
        


}

void Migration::receiver_to_source_propagation()
{



}

void Migration::image_compensation()
{
    modeling->get_runtime();




}

void Migration::export_seismic_volume()
{
    output_image_file = output_image_folder + "image_etc";

    export_binary_float(output_image_file, image, modeling->nPoints);
}