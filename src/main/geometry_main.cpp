# include "../geometry/regular/regular.hpp"
# include "../geometry/circular/circular.hpp"

int main(int argc, char **argv)
{   
    std::vector<Geometry *> geometry = 
    {
        new Regular(), 
        new Circular()
    };

    auto file = std::string(argv[1]);
    auto type = std::stoi(catch_parameter("geometry_type", file));

    geometry[type]->file = file;

    geometry[type]->set_geometry();

    std::vector<Geometry *>().swap(geometry); 

    return 0;
}