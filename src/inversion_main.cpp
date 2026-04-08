# include "../src/inversion/tomography_iso.hpp"
// # include "../src/inversion/tomography_ani.hpp"

int main(int argc, char **argv)
{
    Inversion * inversion = new Tomography_ISO(); 
    
    auto file = std::string(argv[1]);
    auto type = std::stoi(catch_parameter("inv_type", file));

    inversion->parameters = file;

    inversion->set_parameters();
    inversion->import_obsData();

    auto ti = std::chrono::system_clock::now();

    while (true)
    {
        inversion->forward_modeling();
        inversion->check_convergence();

        if (inversion->converged) break; 

        inversion->optimization();
        inversion->model_update();
    }

    auto tf = std::chrono::system_clock::now();

    inversion->export_results();

    std::chrono::duration<double> elapsed_seconds = tf - ti;
    std::cout << "\nRun time: " << elapsed_seconds.count() << " s." << std::endl;

    return 0;
}