# include "../modeling/eikonal/FIM.cuh"
# include "../modeling/eikonal/FSM.cuh"
# include "../modeling/eikonal/PAL.cuh"
# include "../modeling/eikonal/eikonal.hpp"

# include "../modeling/scalar/scalar.hpp"
# include "../modeling/acoustic/acoustic.hpp"
# include "../modeling/elastic/elastic.hpp"

int main(int argc, char **argv)
{
    Modeling * modeling[] = 
    { 
        new Eikonal_pal(), new Eikonal_fim(), new Eikonal_fsm(), 
        new Scalar(), 
        new Acoustic(), 
        new Elastic() 
    }; 

    auto file = std::string(argv[1]);

    auto type = std::stoi(catch_parameter("modeling_type", file)); 

    modeling[type]->set_parameters(file);

    modeling[type]->set_components();

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