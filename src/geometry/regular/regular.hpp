# ifndef REGULAR_HPP
# define REGULAR_HPP

# include "../geometry.hpp"

class Regular : public Geometry
{
    bool reciprocity;

    std::vector<float> shots_SW;
    std::vector<float> shots_NW;
    std::vector<float> shots_SE;

    std::vector<float> nodes_SW;
    std::vector<float> nodes_NW;
    std::vector<float> nodes_SE;

    std::vector<int> shots_nlines;
    std::vector<int> nodes_nlines;
    
    void build_shots();
    void build_nodes();

public:

    void set_parameters(std::string file);
};

# endif