# include "../inversion/tomography/least_squares/least_squares.cuh"
# include "../inversion/tomography/adjoint_state/adjoint_state.cuh"

# include "../inversion/waveform/scalar_isotropic_fwi/scalar_isotropic_fwi.cuh"

int main(int argc, char **argv)
{
    std::vector<Inversion *> inversion = 
    {
        new Least_Squares(), 
        new Adjoint_State(), 
    }; 
    
    auto file = std::string(argv[1]);
    auto type = std::stoi(catch_parameter("inversion_type", file));

    inversion[type]->file = file;

    inversion[type]->set_parameters();
    inversion[type]->import_obs_data();

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