# include "../modeling/eikonal_equation/podvin_and_lecomte/podvin_and_lecomte.cuh"
# include "../modeling/eikonal_equation/fast_sweeping_method/fast_sweeping_method.cuh"
# include "../modeling/eikonal_equation/fast_iterative_method/fast_iterative_method.cuh"

# include "../modeling/wave_equation/scalar/scalar.cuh" 
# include "../modeling/wave_equation/acoustic/acoustic.cuh" 
# include "../modeling/wave_equation/elastic/elastic.cuh" 

int main(int argc, char **argv)
{
    std::vector<Modeling *> modeling = 
    {
        new Podvin_and_Lecomte(),
        new Fast_Iterative_Method(),
        new Fast_Sweeping_Method(),
        
        new Scalar(),
        new Acoustic(),
        new Elastic(), 
    };
    
    auto file = std::string(argv[1]);
    auto type = std::stoi(catch_parameter("modeling_type", file));

    modeling[type]->file = file;

    modeling[type]->set_parameters();
    modeling[type]->set_runtime();

    for (int shot = 0; shot < modeling[type]->total_shots; shot++)
    {
        modeling[type]->shot_id = shot;

        modeling[type]->info_message();
        modeling[type]->initial_setup();
        modeling[type]->forward_solver();
        modeling[type]->build_outputs();
        modeling[type]->export_outputs();
    }

    modeling[type]->get_runtime();
    modeling[type]->free_space();    

    return 0;
}