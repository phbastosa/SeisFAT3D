# include "../inversion/waveform/waveform.hpp"
# include "../inversion/tomography/tomography.hpp"

int main(int argc, char **argv)
{
    Inversion * inversion[] = { new Waveform(), new Tomography() }; 

    auto file = std::string(argv[1]);

    auto type = std::stoi(catch_parameter("inversion_type", file));

    inversion[type]->file = file;
    inversion[type]->set_parameters();

    // inversion[type]->forward_modeling();

    return 0;
}