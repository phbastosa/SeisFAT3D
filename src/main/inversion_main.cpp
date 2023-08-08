# include "../inversion/inversion.hpp"

# include "../inversion/waveform/scalar_fwi/scalar_fwi.hpp"

# include "../inversion/tomography/least_squares/least_squares.hpp"
# include "../inversion/tomography/adjoint_state/adjoint_state.cuh"

int main(int argc, char **argv)
{
    std::vector<Inversion *> inversion = 
    {
        new Least_Squares(), 
        new Adjoint_State(), 
        
        new Scalar_FWI()
    }; 

    // inversion->file = std::string(argv[1]);

    // inversion->set_parameters();
    // inversion->import_obs_data();

    // while (true)
    // {
    //     inversion->forward_modeling();
    //     inversion->compute_residuals();
    //     inversion->check_convergence();

    //     if (inversion->converged) break;

    // //     inversion->optimization();
    // //     inversion->model_update();
    // }

    // inversion->export_results();

    return 0;
}