# include "../modeling/eikonal_equation/isotropic/classical.cuh"
# include "../modeling/eikonal_equation/isotropic/block_FIM.cuh"
# include "../modeling/eikonal_equation/isotropic/ultimate_FSM.cuh"

# include "../modeling/wave_equation/isotropic/elastic.cuh"
# include "../modeling/wave_equation/isotropic/acoustic.cuh"

int main(int argc, char **argv)
{
    std::vector<Modeling *> modeling = 
    {
        new Classical(),
        new Block_FIM(),
        new Ultimate_FSM(),

        new Acoustic(),
        new Elastic()
    };
    
    auto file = std::string(argv[1]);
    auto type = std::stoi(catch_parameter("modeling_type", file));

    modeling[type]->file = file;

    modeling[type]->set_parameters();

    modeling[type]->set_runtime();

    for (int shot = 0; shot < modeling[type]->total_shots; shot++)
    {
        modeling[type]->shot_index = shot;

        modeling[type]->get_information();
        modeling[type]->set_configuration();
        modeling[type]->set_forward_solver();
        
        modeling[type]->export_outputs();
    }

    modeling[type]->get_runtime();
    modeling[type]->free_space();    

    return 0;
}