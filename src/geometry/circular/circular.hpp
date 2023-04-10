# ifndef CIRCULAR_HPP
# define CIRCULAR_HPP

# include "../geometry.hpp"

class Circular : public Geometry
{


public:

    void set_geometry(std::string file, std::string name);
};

# endif