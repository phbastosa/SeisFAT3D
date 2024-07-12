# include "../migration/kirchhoff/kirchhoff.cuh"

int main(int argc, char **argv)
{
    Migration * migration = new Kirchhoff();
    
    auto file = std::string(argv[1]);
    auto type = std::stoi(catch_parameter("modeling_type", file));

    migration->file = file;

    migration->set_parameters();

    migration->image_building();

    return 0;
}