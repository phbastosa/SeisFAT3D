# include "migration/kirchhoff.cuh"

int main(int argc, char **argv)
{
    Migration * migration = new Kirchhoff();
    
    migration->parameters = std::string(argv[1]);

    migration->set_parameters();

    auto ti = std::chrono::system_clock::now();

    migration->image_building();

    auto tf = std::chrono::system_clock::now();

    migration->export_outputs();

    std::chrono::duration<double> elapsed_seconds = tf - ti;
    std::cout << "\nRun time: " << elapsed_seconds.count() << " s." << std::endl;
    
    return 0;
}