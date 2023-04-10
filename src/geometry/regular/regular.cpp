# include "regular.hpp"

void Regular::set_geometry(std::string file, std::string name)
{
    reciprocity = str2bool(catch_parameter("reciprocity", file));    
    import_geometry = str2bool(catch_parameter("import_geometry", file));
    geometry_file = catch_parameter("geometry_folder", file);
 
    splitted = split(catch_parameter(name + "_nlines", file), ',');
    for (auto key : splitted) nlines.push_back(std::stoi(key));

    splitted = split(catch_parameter(name + "_SW", file), ',');
    for (auto key : splitted) SW.push_back(std::stof(key));

    splitted = split(catch_parameter(name + "_NW", file), ',');
    for (auto key : splitted) NW.push_back(std::stof(key));

    splitted = split(catch_parameter(name + "_SE", file), ',');
    for (auto key : splitted) SE.push_back(std::stof(key));

    build_geometry();

    for (int k = 0; k < total; k++)
        std::cout<<z[k]<<" "<<x[k]<<" "<<y[k]<<std::endl;
}

void Regular::build_geometry()
{
    if (import_geometry)
    {
        import_coordinates();
    }
    else
    {
        if (nlines[0] == 1)
        {
            if (nlines[1] == 1)
            {
                if (nlines[2] == 1)
                {
                    total = nlines[0]; // Point source geometry

                    z = new float[total]();
                    x = new float[total]();
                    y = new float[total]();

                    if ((SW[0] == NW[0]) && (SW[0] == SE[0]))
                    {
                        z[0] = SW[0];        
                    }
                    else throw std::invalid_argument("\033[31mGeometry error: wrong z argument for point source selection!\033[0;0m");          
                    
                    if ((SW[1] == NW[1]) && (SW[1] == SE[1]))    
                    {
                        x[0] = SW[1];        
                    }
                    else throw std::invalid_argument("\033[31mGeometry error: wrong x argument for point source selection!\033[0;0m");
                    
                    if ((SW[2] == NW[2]) && (SW[2] == SE[2]))    
                    {
                        y[0] = SW[2];        
                    }
                    else throw std::invalid_argument("\033[31mGeometry error: wrong y argument for point source selection!\033[0;0m");                       
                }
                else
                {
                    total = nlines[2]; // Y line geometry
                    
                    z = new float[total]();
                    x = new float[total]();
                    y = new float[total]();
                
                    if ((SW[0] == NW[0]) && (SW[0] == SE[0]))
                    {
                        for (int k = 0; k < total; k++)
                            z[k] = SW[0];        
                    }
                    else throw std::invalid_argument("\033[31mGeometry error: wrong z argument for point source selection!\033[0;0m");          
                    
                    if ((SW[1] == NW[1]) && (SW[1] == SE[1]))    
                    {
                        for (int k = 0; k < total; k++)
                            x[k] = SW[1];        
                    }
                    else throw std::invalid_argument("\033[31mGeometry error: wrong x argument for point source selection!\033[0;0m");

                    if ((SW[2] == SE[2]) && (SW[2] < NW[2]))
                    {
                        std::vector<float> Y = linspace(SW[2], NW[2], total);            
                        
                        for (int k = 0; k < total; k++)
                            y[k] = Y[k];

                        std::vector<float>().swap(Y);    
                    }
                    else throw std::invalid_argument("\033[31mGeometry error: wrong y argument for point source selection!\033[0;0m");
                }    
            }
            else
            {
                if (nlines[2] == 1)
                {
                    total = nlines[1]; // X line geometry
                    
                    z = new float[total]();
                    x = new float[total]();
                    y = new float[total]();
                
                    if ((SW[0] == NW[0]) && (SW[0] == SE[0]))
                    {
                        for (int k = 0; k < total; k++)
                            z[k] = SW[0];        
                    }
                    else throw std::invalid_argument("\033[31mGeometry error: wrong z argument for point source selection!\033[0;0m");          
                    
                    if ((SW[1] == NW[1]) && (SW[1] < SE[1]))    
                    {
                        std::vector<float> X = linspace(SW[1], SE[1], total);            

                        for (int k = 0; k < total; k++)
                            x[k] = X[k];        

                        std::vector<float>().swap(X);    
                    }
                    else throw std::invalid_argument("\033[31mGeometry error: wrong x argument for point source selection!\033[0;0m");

                    if ((SW[2] == NW[2]) && (SW[2] == SE[2]))
                    {
                        for (int k = 0; k < total; k++)
                            y[k] = SW[2];
                    }
                    else throw std::invalid_argument("\033[31mGeometry error: wrong y argument for point source selection!\033[0;0m");                
                }
                else
                {
                    std::cout << "XY plane geometry" << std::endl;
                
                
                
                
                }    
            }
        }
        else
        {
            if (nlines[1] == 1)
            {
                if (nlines[2] == 1)
                {
                    std::cout << "Z line geometry" << std::endl;
                    total = nlines[0]; // Z line geometry
                    
                    z = new float[total]();
                    x = new float[total]();
                    y = new float[total]();
                
                    if ((SW[0] == SE[0]) && (SW[0] < NW[0]))
                    {
                        std::vector<float> Z = linspace(SW[0], NW[0], total);            

                        for (int k = 0; k < total; k++)
                            z[k] = Z[k];        

                        std::vector<float>().swap(Z);    
                    }
                    else throw std::invalid_argument("\033[31mGeometry error: wrong z argument for point source selection!\033[0;0m");          
                    
                    if ((SW[1] == NW[1]) && (SW[1] == SE[1]))    
                    {
                        for (int k = 0; k < total; k++)
                            x[k] = SW[1];        
                    }
                    else throw std::invalid_argument("\033[31mGeometry error: wrong x argument for point source selection!\033[0;0m");

                    if ((SW[2] == NW[2]) && (SW[2] == SE[2]))
                    {
                        for (int k = 0; k < total; k++)
                            y[k] = SW[2];
                    }
                    else throw std::invalid_argument("\033[31mGeometry error: wrong y argument for point source selection!\033[0;0m");                                
                }
                else
                {
                    std::cout << "ZY plane geometry" << std::endl;
                
                
                
                
                }    
            }
            else
            {
                if (nlines[2] == 1)
                {
                    std::cout << "ZX plane geometry" << std::endl;
                
                
                
                }
                else
                {
                    std::cout << "ZXY cube geometry" << std::endl;
                
                
                
                }    
            }
        }
    }
}

std::vector<float> Regular::linspace(float xi, float xf, int n)
{
    std::vector<float> linspaced;
    
    if (n == 0) return linspaced;
    if (n == 1)
    {
        linspaced.push_back(xi);
        return linspaced;
    } 

    linspaced.reserve(n);

    float delta = (xf - xi) / (n - 1);

    for (int i = 0; i < n; i++)
    {
        linspaced.emplace_back(xi + (float)(delta*i));
    }

    return linspaced;
}

