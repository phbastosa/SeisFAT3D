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


    return 0;
}