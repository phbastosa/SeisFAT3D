# ifndef REGULAR_HPP
# define REGULAR_HPP

# include "../geometry.hpp"

class Regular : public Geometry
{
    void set_name();

    void set_parameters(std::string file);
};

# endif