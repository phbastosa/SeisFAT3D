# include "geometry/regular/regular.hpp"
# include "geometry/circular/circular.hpp"

int main(int argc, char **argv)
{
    Geometry * geometry[] = { new Regular(), new Circular()}; 

    std::string file = std::string(argv[1]);

    int type = std::stoi(catch_parameter("geometry_type", file));

    geometry[type]->set_geometry(file);

    return 0;
}