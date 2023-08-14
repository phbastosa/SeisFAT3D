# include "geometry.hpp"

void Geometry::set_general_parameters()
{
    reciprocity = str2bool(catch_parameter("reciprocity", file));
    import_geometry = str2bool(catch_parameter("import_geometry", file));

    shots_file = catch_parameter("shots_file", file);
    nodes_file = catch_parameter("nodes_file", file);
}

void Geometry::set_reciprocity()
{
    std::swap(shots.x, nodes.x);    
    std::swap(shots.y, nodes.y);    
    std::swap(shots.z, nodes.z);

    std::swap(shots.total, nodes.total);
}

void Geometry::import_coordinates()
{
    std::vector<std::string> elements;

    import_text_file(shots_file, elements);

    shots.total = elements.size();

    shots.x = new float[shots.total]();
    shots.y = new float[shots.total]();
    shots.z = new float[shots.total]();

    for (int i = 0; i < shots.total; i++)
    {
        splitted = split(elements[i], ',');

        shots.x[i] = std::stof(splitted[0]);
        shots.y[i] = std::stof(splitted[1]);
        shots.z[i] = std::stof(splitted[2]);
    }    

    std::vector<std::string>().swap(elements);

    import_text_file(nodes_file, elements);

    nodes.total = elements.size();

    nodes.x = new float[nodes.total]();
    nodes.y = new float[nodes.total]();
    nodes.z = new float[nodes.total]();

    for (int i = 0; i < nodes.total; i++)
    {
        splitted = split(elements[i], ',');

        nodes.x[i] = std::stof(splitted[0]);
        nodes.y[i] = std::stof(splitted[1]);
        nodes.z[i] = std::stof(splitted[2]);
    }    

    std::vector<std::string>().swap(elements);
}

void Geometry::export_coordinates()
{
    auto folder = std::string("../inputs/geometry/");
    auto shots_path = folder + std::string("xyz_shot_positions.txt");
    auto nodes_path = folder + std::string("xyz_node_positions.txt");

    std::ofstream sfile(shots_path, std::ios::out);        

    if (sfile.is_open()) 
    {    
        for (int i = 0; i < shots.total; i++)        
            sfile <<shots.x[i]<<", "<<shots.y[i]<<", "<<shots.z[i]<<std::endl;    
        
    }
    else
    {
        throw std::invalid_argument("Error: file " + shots_path + " could not be opened!");
    }

    std::cout<<"File " + shots_path + " was succesfully written."<<std::endl;

    sfile.close();

    std::ofstream nfile(nodes_path, std::ios::out);        

    if (nfile.is_open()) 
    {    
        for (int i = 0; i < nodes.total; i++)        
        {   
            nfile <<nodes.x[i]<<", "<<nodes.y[i]<<", "<<nodes.z[i]<<std::endl;    
        }
    }
    else
    {
        throw std::invalid_argument("Error: file " + nodes_path + " could not be opened!");
    }

    std::cout<<"File " + nodes_path + " was succesfully written."<<std::endl;

    nfile.close();
}

void Geometry::set_regular(Coord &obj)
{
    if (nlines[0] == 1)
    {
        if (nlines[1] == 1)
        {
            if (nlines[2] == 1)
            {
                obj.total = nlines[0]; // Point source geometry

                obj.z = new float[obj.total]();
                obj.x = new float[obj.total]();
                obj.y = new float[obj.total]();

                if ((SW[0] == NW[0]) && (SW[0] == SE[0]))
                {
                    obj.z[0] = SW[0];        
                }
                else throw std::invalid_argument("\033[31mGeometry error: wrong z argument for point selection!\033[0;0m");          
                
                if ((SW[1] == NW[1]) && (SW[1] == SE[1]))    
                {
                    obj.x[0] = SW[1];        
                }
                else throw std::invalid_argument("\033[31mGeometry error: wrong x argument for point selection!\033[0;0m");
                
                if ((SW[2] == NW[2]) && (SW[2] == SE[2]))    
                {
                    obj.y[0] = SW[2];        
                }
                else throw std::invalid_argument("\033[31mGeometry error: wrong y argument for point selection!\033[0;0m");                       
            }
            else
            {
                obj.total = nlines[2]; // Y line geometry
                
                obj.z = new float[obj.total]();
                obj.x = new float[obj.total]();
                obj.y = new float[obj.total]();
            
                if ((SW[0] == NW[0]) && (SW[0] == SE[0]))
                {
                    for (int i = 0; i < obj.total; i++)
                        obj.z[i] = SW[0];        
                }
                else throw std::invalid_argument("\033[31mGeometry error: wrong z argument for line selection!\033[0;0m");          
                
                if ((SW[1] == NW[1]) && (SW[1] == SE[1]))    
                {
                    for (int j = 0; j < obj.total; j++)
                        obj.x[j] = SW[1];        
                }
                else throw std::invalid_argument("\033[31mGeometry error: wrong x argument for line selection!\033[0;0m");

                if ((SW[2] == SE[2]) && (SW[2] < NW[2]))
                {
                    std::vector<float> Y = linspace(SW[2], NW[2], obj.total);            
                    
                    for (int k = 0; k < obj.total; k++)
                        obj.y[k] = Y[k];

                    std::vector<float>().swap(Y);    
                }
                else throw std::invalid_argument("\033[31mGeometry error: wrong y argument for line selection!\033[0;0m");
            }    
        }
        else
        {
            if (nlines[2] == 1)
            {
                obj.total = nlines[1]; // X line geometry
                
                obj.z = new float[obj.total]();
                obj.x = new float[obj.total]();
                obj.y = new float[obj.total]();
            
                if ((SW[0] == NW[0]) && (SW[0] == SE[0]))
                {
                    for (int i = 0; i < obj.total; i++)
                        obj.z[i] = SW[0];        
                }
                else throw std::invalid_argument("\033[31mGeometry error: wrong z argument for line selection!\033[0;0m");          
                
                if ((SW[1] == NW[1]) && (SW[1] < SE[1]))    
                {
                    std::vector<float> X = linspace(SW[1], SE[1], obj.total);            

                    for (int j = 0; j < obj.total; j++)
                        obj.x[j] = X[j];        

                    std::vector<float>().swap(X);    
                }
                else throw std::invalid_argument("\033[31mGeometry error: wrong x argument for line selection!\033[0;0m");

                if ((SW[2] == NW[2]) && (SW[2] == SE[2]))
                {
                    for (int k = 0; k < obj.total; k++)
                        obj.y[k] = SW[2];
                }
                else throw std::invalid_argument("\033[31mGeometry error: wrong y argument for line selection!\033[0;0m");                
            }
            else
            {
                obj.total = nlines[1] * nlines[2]; // XY plane geometry
                
                obj.z = new float[obj.total]();
                obj.x = new float[obj.total]();
                obj.y = new float[obj.total]();
            
                if ((SW[0] == NW[0]) && (SW[0] == SE[0]))
                {
                    for (int i = 0; i < obj.total; i++)
                        obj.z[i] = SW[0];        
                }
                else throw std::invalid_argument("\033[31mGeometry error: wrong z argument plane selection!\033[0;0m");          
                
                if ((SW[1] == NW[1]) && (SW[1] < SE[1]))    
                {
                    if ((SW[2] == SE[2]) && (SW[2] < NW[2]))    
                    {
                        std::vector<float> X = linspace(SW[1], SE[1], nlines[1]);            
                        std::vector<float> Y = linspace(SW[2], NW[2], nlines[2]);            

                        for (int k = 0; k < Y.size(); k++)
                        {
                            for (int j = 0; j < X.size(); j++)
                            {
                                obj.x[k + j*nlines[2]] = X[j];
                                obj.y[k + j*nlines[2]] = Y[k];
                            }
                        }  
                    
                        std::vector<float>().swap(X);    
                        std::vector<float>().swap(Y);    
                    }        
                }
                else throw std::invalid_argument("\033[31mGeometry error: wrong xy argument for plane selection!\033[0;0m");
            }    
        }
    }
    else
    {
        if (nlines[1] == 1)
        {
            if (nlines[2] == 1)
            {
                obj.total = nlines[0]; // Z line geometry
                
                obj.z = new float[obj.total]();
                obj.x = new float[obj.total]();
                obj.y = new float[obj.total]();
            
                if ((SW[0] == SE[0]) && (SW[0] < NW[0]))
                {
                    std::vector<float> Z = linspace(SW[0], NW[0], obj.total);            

                    for (int i = 0; i < obj.total; i++)
                        obj.z[i] = Z[i];        

                    std::vector<float>().swap(Z);    
                }
                else throw std::invalid_argument("\033[31mGeometry error: wrong z argument for line selection!\033[0;0m");          
                
                if ((SW[1] == NW[1]) && (SW[1] == SE[1]))    
                {
                    for (int j = 0; j < obj.total; j++)
                        obj.x[j] = SW[1];        
                }
                else throw std::invalid_argument("\033[31mGeometry error: wrong x argument for line selection!\033[0;0m");

                if ((SW[2] == NW[2]) && (SW[2] == SE[2]))
                {
                    for (int k = 0; k < obj.total; k++)
                        obj.y[k] = SW[2];
                }
                else throw std::invalid_argument("\033[31mGeometry error: wrong y argument for line selection!\033[0;0m");                                
            }
            else
            {
                obj.total = nlines[0] * nlines[2]; // ZY plane geometry
                
                obj.z = new float[obj.total]();
                obj.x = new float[obj.total]();
                obj.y = new float[obj.total]();
            
                if ((SW[1] == NW[1]) && (SW[1] == SE[1]))
                {
                    for (int j = 0; j < obj.total; j++)
                        obj.x[j] = SW[0];        
                }
                else throw std::invalid_argument("\033[31mGeometry error: wrong x argument plane selection!\033[0;0m");          
                
                if ((SW[0] == SE[0]) && (SW[0] < NW[0]))    
                {
                    if ((SW[2] == SE[2]) && (SW[2] < NW[2]))    
                    {
                        std::vector<float> Z = linspace(SW[0], NW[0], nlines[0]);            
                        std::vector<float> Y = linspace(SW[2], NW[2], nlines[2]);            

                        for (int i = 0; i < Z.size(); i++)
                        {
                            for (int k = 0; k < Y.size(); k++)
                            {
                                obj.z[i + k*nlines[0]] = Z[i];
                                obj.y[i + k*nlines[0]] = Y[k];
                            }
                        }  
                    
                        std::vector<float>().swap(Z);    
                        std::vector<float>().swap(Y);    
                    }        
                }
                else throw std::invalid_argument("\033[31mGeometry error: wrong zy argument for plane selection!\033[0;0m");                
            }    
        }
        else
        {
            if (nlines[2] == 1)
            {
                obj.total = nlines[0] * nlines[1]; // ZX plane geometry
                
                obj.z = new float[obj.total]();
                obj.x = new float[obj.total]();
                obj.y = new float[obj.total]();
            
                if ((SW[2] == NW[2]) && (SW[2] == SE[2]))
                {
                    for (int k = 0; k < obj.total; k++)
                        obj.y[k] = SW[2];        
                }
                else throw std::invalid_argument("\033[31mGeometry error: wrong y argument plane selection!\033[0;0m");          
                
                if ((SW[0] == SE[0]) && (SW[0] < NW[0]))    
                {
                    if ((SW[1] == NW[1]) && (SW[1] < SE[1]))    
                    {
                        std::vector<float> Z = linspace(SW[0], NW[0], nlines[0]);            
                        std::vector<float> X = linspace(SW[1], SE[1], nlines[1]);            

                        for (int i = 0; i < Z.size(); i++)
                        {
                            for (int j = 0; j < X.size(); j++)
                            {
                                obj.z[i + j*nlines[0]] = Z[i];
                                obj.x[i + j*nlines[0]] = X[j];
                            }
                        }  
                    
                        std::vector<float>().swap(Z);    
                        std::vector<float>().swap(X);    
                    }        
                }
                else throw std::invalid_argument("\033[31mGeometry error: wrong zy argument for plane selection!\033[0;0m");                                
            }
            else throw std::invalid_argument("\033[31mGeometry error: xyz cube selection is not available!\033[0;0m");
        }
    }

    std::vector<int>().swap(nlines);
  
    std::vector<float>().swap(NW);
    std::vector<float>().swap(SW);
    std::vector<float>().swap(SE);
}

std::vector<float> Geometry::linspace(float xi, float xf, int n)
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
