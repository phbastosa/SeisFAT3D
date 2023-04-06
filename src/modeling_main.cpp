# include "modeling/eikonal/eikonal.hpp"
# include "modeling/scalar/scalar.hpp"
# include "modeling/acoustic/acoustic.hpp"
# include "modeling/elastic/elastic.hpp"

int main(int argc, char **argv)
{
    Modeling * modeling[] =
    {
        new Eikonal(),
        new Scalar(),
        new Acoustic(),
        new Elastic()
    }; 

    for (int type = 0; type < 4; type++)
    {
        modeling[type]->set_name();
        modeling[type]->get_name();
    }

    return 0;
}