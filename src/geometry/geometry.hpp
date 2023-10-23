# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include "../io/io.hpp"

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
    void set_general_parameters();

    void set_regular(Coord &obj);

public:

    Coord shots;
    Coord nodes;

    std::string file;

    virtual void set_geometry() = 0;     
};

# endif
