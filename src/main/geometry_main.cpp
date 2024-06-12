# include "../geometry/geometry.hpp"

int main(int argc, char **argv)
{   
    auto * geometry = new Geometry();

    geometry->file = std::string(argv[1]);

    geometry->set_geometry();

    return 0;
}