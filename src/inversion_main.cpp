# include "inversion/waveform/waveform.hpp"
# include "inversion/tomography/tomography.hpp"

int main(int argc, char **argv)
{
    Inversion * inversion[] =
    {
        new Waveform(),
        new Tomography()
    }; 

    for (int type = 0; type < 2; type++)
    {
        inversion[type]->set_name();
        inversion[type]->get_name();
    }

    return 0;
}