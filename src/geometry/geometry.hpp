# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include "../utils/input_output/io.hpp"

class Geometry
{
protected:

    std::string name;        
    std::string geometry_file;

public:
    
    void get_name();
    virtual void set_name() = 0;

    float * x = nullptr;
    float * y = nullptr;
    float * z = nullptr;

    int * beg_relative = nullptr;
    int * end_relative = nullptr;

    void import_coordinates();
    void export_coordinates();

    virtual void set_parameters(std::string file) = 0;     
};

# endif