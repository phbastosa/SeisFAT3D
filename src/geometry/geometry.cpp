# include "geometry.hpp"

void Geometry::set_parameters()
{    
    std::vector<std::string> SPS;
    std::vector<std::string> RPS;
    std::vector<std::string> AUX;

    std::string SPS_file = catch_parameter("SPS", parameters);
    std::string RPS_file = catch_parameter("RPS", parameters);

    import_text_file(SPS_file, SPS); 
    import_text_file(RPS_file, RPS); 

    nsrc = SPS.size();
    nrec = RPS.size();

    xsrc = new float[nsrc]();
    ysrc = new float[nsrc]();
    zsrc = new float[nsrc]();

    xrec = new float[nrec]();
    yrec = new float[nrec]();
    zrec = new float[nrec]();

    for (int i = 0; i < nsrc; i++)
    {
        AUX = split(SPS[i], ',');

        xsrc[i] = std::stof(AUX[0]);
        ysrc[i] = std::stof(AUX[1]);
        zsrc[i] = std::stof(AUX[2]);
    }    

    std::vector<std::string>().swap(AUX);

    for (int i = 0; i < nrec; i++)
    {
        AUX = split(RPS[i], ',');

        xrec[i] = std::stof(AUX[0]);
        yrec[i] = std::stof(AUX[1]);
        zrec[i] = std::stof(AUX[2]);
    }    
    
    std::vector<std::string>().swap(AUX);
}