# include "geometry/regular/regular.hpp"
# include "geometry/circular/circular.hpp"
# include "geometry/streamer/streamer.hpp"

int main(int argc, char **argv)
{
    Geometry * shots[] = { new Regular(), new Circular(), new Streamer() }; 
    Geometry * nodes[] = { new Regular(), new Circular(), new Streamer() }; 

    std::string file = std::string(argv[1]);

    int type = std::stoi(catch_parameter("geometry_type", file));

    switch (type)
    {
        case 0:
            shots[0]->set_geometry(file, "shots");
            nodes[0]->set_geometry(file, "nodes");
            break;

        case 1:
            shots[1]->set_geometry(file, "shots");
            nodes[0]->set_geometry(file, "nodes");
            break;

        case 2:
            shots[0]->set_geometry(file, "shots");
            nodes[2]->set_geometry(file, "nodes");
            break;

        default:
            throw std::invalid_argument("\033[31mGeometry error: wrong type!\033[0;0m");          
            break;
    }

    return 0;
}