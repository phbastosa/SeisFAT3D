# include "eikonal.hpp"

void Eikonal::set_parameters(std::string file)
{
    Modeling::set_parameters(file);

    vp_model_file = catch_parameter("vp_model_file", file);


}