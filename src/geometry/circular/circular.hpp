# ifndef CIRCULAR_HPP
# define CIRCULAR_HPP

# include "../geometry.hpp"

class Circular : public Geometry
{
private:
    
    int spacing;
    float xc, yc, zc;

    std::vector<float> offsets;

    void set_circular();

public:

    void set_geometry();
};

# endif