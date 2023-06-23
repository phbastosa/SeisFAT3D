# include "../modeling/modeling.cuh"

int main(int argc, char **argv)
{
    Modeling * modeling = new Modeling();

    auto file = std::string(argv[1]);

    modeling->file = file;

    modeling->set_parameters();
    modeling->set_slowness();
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

    return 0;
}