# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include "../utils/input_output/io.hpp"

class Coord
{
public:

    int total;

    float * x = nullptr;
    float * y = nullptr;
    float * z = nullptr;    
};

class Geometry
{
private:

    std::vector<float> linspace(float xi, float xf, int n);

protected:

    bool reciprocity;
    
    bool import_geometry;
    std::string shots_file;
    std::string nodes_file;

    std::vector<int> nlines;
    std::vector<float> SW, NW, SE;    
    std::vector<std::string> splitted;

    void set_reciprocity();
    void import_coordinates();
    void export_coordinates();

    void set_regular(Coord &obj);

public:

    Coord shots;
    Coord nodes;

    virtual void set_geometry(std::string file);     
};

# endif
