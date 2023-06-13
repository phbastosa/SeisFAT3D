# include "tomography.hpp"

void Tomography::set_parameters()
{
    modeling = new Eikonal();

    modeling->file = file;

    modeling->set_parameters();

    

}

void Tomography::import_obs_data()
{


}

void Tomography::forward_modeling()
{
    modeling->set_components();

    modeling->set_runtime();

    for (int shot = 0; shot < modeling->total_shots; shot++)
    {
        modeling->shot_id = shot;

        modeling->info_message();

        modeling->initial_setup();

        modeling->forward_solver();

        modeling->build_outputs();

        modeling->export_outputs();
    }

    modeling->get_runtime();

    modeling->free_space(); 
}

void Tomography::compute_residuals()
{



}

void Tomography::check_convergence()
{



}

void Tomography::optimization()
{


}

void Tomography::model_update()
{



}

