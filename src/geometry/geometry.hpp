# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include "../utils/input_output/io.hpp"

class Geometry
{
protected:

    std::string name;        

    bool import_geometry;
    std::string geometry_folder;

public:
    
    float * shot_x = nullptr;
    float * shot_y = nullptr;
    float * shot_z = nullptr;

    float * node_x = nullptr;
    float * node_y = nullptr;
    float * node_z = nullptr;

    int * beg_relative = nullptr;
    int * end_relative = nullptr;

    void import_coordinates();
    void export_coordinates();

    virtual void set_parameters(std::string file) = 0;     
};

# endif