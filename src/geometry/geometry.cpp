# include "geometry.hpp"

void Geometry::set_parameters()
{    
    std::vector<std::string> SPS;
    std::vector<std::string> RPS;
    std::vector<std::string> XPS;
    std::vector<std::string> AUX;

    std::string SPS_file = catch_parameter("SPS", parameters);
    std::string RPS_file = catch_parameter("RPS", parameters);
    std::string XPS_file = catch_parameter("XPS", parameters);

    import_text_file(SPS_file, SPS); 
    import_text_file(RPS_file, RPS); 
    import_text_file(XPS_file, XPS); 

    nsrc = SPS.size();
    nrec = RPS.size();
    nrel = XPS.size();

    sInd = new int[nrel]();
    iRec = new int[nrel]();
    fRec = new int[nrel]();

    spread = new int[nrel]();

    xsrc = new float[nsrc]();
    ysrc = new float[nsrc]();
    zsrc = new float[nsrc]();

    xrec = new float[nrec]();
    yrec = new float[nrec]();
    zrec = new float[nrec]();

    for (int i = 0; i < nrel; i++)
    {
        AUX = split(XPS[i], ',');

        sInd[i] = std::stoi(AUX[0]);
        iRec[i] = std::stoi(AUX[1]);
        fRec[i] = std::stoi(AUX[2]);

        spread[i] = fRec[i] - iRec[i];
    }    

    std::vector<std::string>().swap(AUX);

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
}