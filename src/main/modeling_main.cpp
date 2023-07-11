# include "../modeling/modeling.hpp"

# include "../modeling/high_frequency/eikonal/podvin_and_lecomte/podvin_and_lecomte.cuh"
# include "../modeling/high_frequency/eikonal/fast_sweeping_method/fast_sweeping_method.cuh"
# include "../modeling/high_frequency/eikonal/fast_iterative_method/fast_iterative_method.cuh"

# include "../modeling/low_frequency/scalar/fdm_isotropic/scalar_isotropic.cuh" 
# include "../modeling/low_frequency/acoustic/fdm_isotropic/acoustic_isotropic.cuh" 
# include "../modeling/low_frequency/elastic/fdm_isotropic/elastic_isotropic.cuh" 

int main(int argc, char **argv)
{
    std::vector<Modeling *> modeling = 
    {
        new Podvin_and_Lecomte(),
        new Fast_Iterative_Method(),
        new Fast_Sweeping_Method(),
        
        new Scalar_Isotropic(),
        new Acoustic_Isotropic(),
        new Elastic_Isotropic() 
    };
    
    if (!fileExists(std::string(argv[1])))
        throw std::invalid_argument("\033[31mError: " + std::string(argv[1]) + " could not be opened!\033[0;0m");

    auto file = std::string(argv[1]);

    if (!isInteger(catch_parameter("modeling_type", file)))
        throw std::invalid_argument("\033[31mError: Wrong modeling type! \033[0;0m");
    
    auto type = std::stoi(catch_parameter("modeling_type", file));

    if ((type < 0) || (type >= modeling.size()))
        throw std::invalid_argument("\033[31mError: Modeling type is out of bounds! \033[0;0m");

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