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

    for (int type = 0; type < 3; type++)
    {
        geometry[type]->set_name();
        geometry[type]->get_name();
    }

    return 0;
}