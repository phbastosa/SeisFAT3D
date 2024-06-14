# include "../modeling/hfreq_modeling.cuh"
# include "../modeling/lfreq_modeling.cuh"

int main(int argc, char **argv)
{
    std::vector<Modeling *> modeling = 
    {
        new hfreq_Modeling(), 
        new lfreq_Modeling()
    };
    
    auto file = std::string(argv[1]);
    auto type = std::stoi(catch_parameter("modeling_type", file));

    modeling[type]->file = file;

    modeling[type]->set_parameters();

    modeling[type]->set_runtime();

    for (int shot = 0; shot < modeling[type]->total_shots; shot++)
    {
        modeling[type]->shot_index = shot;

        modeling[type]->print_information();
        
        modeling[type]->set_initial_conditions();

        modeling[type]->forward_propagation();
        
        modeling[type]->export_outputs();
    }

    modeling[type]->get_runtime();

    modeling[type]->free_space();    

    return 0;
}