# ifndef REGULAR_HPP
# define REGULAR_HPP

# include "../geometry.hpp"

class Regular : public Geometry
{
private:
    
    void set_relational();

public:

    void set_geometry(std::string file);
};

# endif