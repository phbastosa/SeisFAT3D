# include "migration/kirchhoff_iso.cuh"
# include "migration/kirchhoff_ani.cuh"

int main(int argc, char **argv)
{
    std::vector<Migration *> migration = 
    {
        new Kirchhoff_ISO(),
        new Kirchhoff_ANI()
    };

    auto file = std::string(argv[1]);
    auto type = std::stoi(catch_parameter("migration_type", file));    

    migration[type]->parameters = file;

    migration[type]->set_parameters();

    auto ti = std::chrono::system_clock::now();

    migration[type]->image_building();

    auto tf = std::chrono::system_clock::now();

    migration[type]->export_outputs();

    std::chrono::duration<double> elapsed_seconds = tf - ti;
    std::cout << "\nRun time: " << elapsed_seconds.count() << " s." << std::endl;
    
    return 0;
}