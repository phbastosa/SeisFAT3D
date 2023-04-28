# include "migration/kirchhoff/kirchhoff.hpp"
# include "migration/reverseTime/reverseTime.hpp"

int main(int argc, char **argv)
{
    Migration * migration[] =
    {
        new Kirchhoff(),
        new ReverseTime(),
    }; 

    for (int type = 0; type < 2; type++)
    {
        migration[type]->set_name();
        migration[type]->get_name();
    }

    return 0;
}