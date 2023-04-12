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
    std::string relational;

    std::vector<int> nlines;
    std::vector<float> SW, NW, SE;    
    std::vector<std::string> splitted;

    void set_reciprocity();
    void import_coordinates();
    void export_coordinates();

    void set_regular(Coord &obj);

    virtual void set_relational() = 0;

public:

    Coord shots;
    Coord nodes;

    int * beg_relation = nullptr;
    int * end_relation = nullptr;

    virtual void set_geometry(std::string file) = 0;     
};

# endif
