# include "geometry/regular/regular.hpp"
# include "geometry/circular/circular.hpp"
# include "geometry/streamer/streamer.hpp"

int main(int argc, char **argv)
{
    Geometry * geometry[] =
    {
        new Regular(),
        new Circular(),
        new Streamer()
    }; 

    std::string file = std::string(argv[1]);

    int type = std::stoi(catch_parameter("geometry_type", file));

    geometry[type]->set_parameters(file);

    // geometry[0]->set_name();
    // geometry[0]->get_name();

    return 0;
}