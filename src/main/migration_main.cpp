# include "../migration/kirchhoff/kirchhoff.cuh"

int main(int argc, char **argv)
{
    Migration * migration = new Kirchhoff();
    
    migration->get_name();

    return 0;
}