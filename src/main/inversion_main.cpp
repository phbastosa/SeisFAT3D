# include "../inversion/tomography/tomography.cuh"

int main(int argc, char **argv)
{
    Inversion * inversion = new Tomography(); 

    inversion->file = std::string(argv[1]);

    inversion->set_parameters();
    inversion->import_obs_data();

    while (true)
    {
        inversion->forward_modeling();
        inversion->compute_residuals();
        inversion->check_convergence();

        if (inversion->converged) break;

    //     inversion->optimization();
    //     inversion->model_update();
    }

    inversion->export_results();

    return 0;
}