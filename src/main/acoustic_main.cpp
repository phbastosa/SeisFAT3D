# include "../seismogram/acoustic.cuh"

int main(int argc, char **argv)
{   
    auto * modeling = new Acoustic();

    modeling->file = std::string(argv[1]);

    modeling->set_patameters();

    modeling->set_runtime();

    for (int shot = 0; shot < modeling->total_shots; shot++)
    {
        modeling->shot_id = shot;
        
        modeling->initial_setup();
        modeling->forward_solver();
        modeling->export_outputs();
    }

    modeling->get_runtime();

    modeling->free_space();

    return 0;
}