# include "io.hpp"

bool str2bool(std::string s)
{
    bool b;

    std::for_each(s.begin(), s.end(), [](char & c){c = ::tolower(c);});
    std::istringstream(s) >> std::boolalpha >> b;

    return b;
}

bool isInteger(const std::string& input) 
{
    std::stringstream ss(input);
    int value;
    return (ss >> value) && (ss.eof());
}

bool fileExists(const std::string& filename) 
{
    std::ifstream file(filename);
    return file.good();
}

void import_binary_float(std::string path, float * array, int n)
{
    std::ifstream file(path, std::ios::in);

    if (file.is_open()) 
    {    
        file.read((char *) array, n * sizeof(float));
    }
    else
    {
        throw std::invalid_argument("\033[31mError: " + path + " could not be opened!\033[0;0m");
    }

    file.close();    
}

void export_binary_float(std::string path, float * array, int n)
{
    std::ofstream file(path, std::ios::out);
    
    if (file.is_open()) 
    {    
        file.write((char *) array, n * sizeof(float));
    }
    else
    {
        throw std::invalid_argument("\033[31mError: " + path + " could not be opened!\033[0;0m");
    }

    std::cout<<"Binary file " + path + " was succesfully written."<<std::endl;

    file.close();
}

void import_text_file(std::string path, std::vector<std::string> &elements)
{
    std::ifstream file(path, std::ios::in);
    
    if (file.is_open()) 
    {    
        std::string line;

        while(getline(file, line))
        {
            elements.push_back(line);
        }
    }
    else
    {
        throw std::invalid_argument("\033[31mError: " + path + " could not be opened!\033[0;0m");
    }

    file.close();
}

void export_text_file(std::string path, std::vector<std::string> &elements)
{
    std::ofstream file(path, std::ios::out);
    
    if (file.is_open()) 
    {    
        for (auto line : elements)
        {
            file << line <<"\n";
        }
    }
    else
    {
        throw std::invalid_argument("\033[31mError: " + path + " could not be opened!\033[0;0m");
    }

    std::cout<<"Text file " + path + " was succesfully written."<<std::endl;

    file.close();
}

std::string catch_parameter(std::string target, std::string file)
{
    char spaces = ' ';
    char comment = '#';

    std::string line;
    std::string variable;

    std::ifstream parameters(file);

    if (parameters.is_open())
    {
        while (getline(parameters, line))
        {           
            if ((line.front() != comment) && (line.front() != spaces))        
            {
                if (line.find(target) == 0)
                {
                    for (int i = line.find("=")+2; i < line.size(); i++)
                    {    
                        if (line[i] == '#') break;
                        variable += line[i];            
                    }

                    break;
                }
            }                 
        }
        parameters.close();
    }        

    // Quality control for file paths

    if (variable.find('"') == 0)
    {
        remove(variable.begin(), variable.end(), '"');
    }
    else if (variable.find("[") == 0)
    {
        remove(variable.begin(), variable.end(), '[');
        remove(variable.begin(), variable.end(), ']');
    }

    variable.erase(remove(variable.begin(), variable.end(), ' '), variable.end());

    return variable;
}

std::vector<std::string> split(std::string s, char delimiter)
{
    std::string token;
    std::vector<std::string> tokens;
    std::istringstream tokenStream(s);

    while (getline(tokenStream, token, delimiter)) 
        tokens.push_back(token);
   
    return tokens;
}