# include "modeling/eikonal/eikonal.cuh"
# include "modeling/scalar/scalar.hpp"
# include "modeling/acoustic/acoustic.hpp"
# include "modeling/elastic/elastic.hpp"

int main(int argc, char **argv)
{
    Modeling * modeling[] =
    {
        new Eikonal(),
        new Scalar(),
        new Acoustic(),
        new Elastic()
    }; 

    int type = std::stoi(catch_parameter("modeling_type", std::string(argv[1]))); 

    modeling[type]->set_parameters(std::string(argv[1]));

    modeling[type]->set_runtime();

    for (int shot = 0; shot < modeling[type]->total_shots; shot++)
    {
        modeling[type]->shot_id = shot;

        


        modeling[type]->export_outputs();
    }

    modeling[type]->get_runtime();

    return 0;
}