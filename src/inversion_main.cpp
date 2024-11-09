# include "inversion/least_squares.hpp"
# include "inversion/adjoint_state.cuh"

int main(int argc, char **argv)
{
    std::vector<Tomography *> inversion = 
    {
        new Least_Squares(), 
        new Adjoint_State(), 
    }; 
    
    auto file = std::string(argv[1]);
    auto type = std::stoi(catch_parameter("inversion_type", file));

    inversion[type]->parameters = file;

    inversion[type]->set_parameters();
    inversion[type]->import_obsData();

    while (true)
    {
        inversion[type]->forward_modeling();
        inversion[type]->check_convergence();

        if (inversion[type]->converged) break; 

        inversion[type]->optimization();
        inversion[type]->model_update();
    }

    inversion[type]->export_results();

    return 0;
}