# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include "../utils/input_output/io.hpp"

class Geometry
{
protected:

    bool import_geometry;
    std::string geometry_file;

public:
    
    int total;

    float * x = nullptr;
    float * y = nullptr;
    float * z = nullptr;

    int * beg_relative = nullptr;
    int * end_relative = nullptr;

    void import_coordinates();
    void export_coordinates();

    virtual void set_geometry(std::string file, std::string name) = 0;     
};

# endif