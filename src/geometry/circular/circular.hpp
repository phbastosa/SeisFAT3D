# ifndef CIRCULAR_HPP
# define CIRCULAR_HPP

# include "../geometry.hpp"

class Circular : public Geometry
{
    void set_name();

    void set_parameters(std::string file);
};

# endif