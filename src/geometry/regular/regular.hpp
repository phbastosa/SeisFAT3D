# ifndef REGULAR_HPP
# define REGULAR_HPP

# include "../geometry.hpp"

class Regular : public Geometry
{
    bool reciprocity;

    std::string parameters;

    std::vector<float> SW;
    std::vector<float> NW;
    std::vector<float> SE;

    std::vector<int> nlines;

    std::vector<std::string> splitted;

    std::vector<float> linspace(float xi, float xf, int n);

    void build_geometry();

public:

    void set_geometry(std::string file, std::string name);
};

# endif