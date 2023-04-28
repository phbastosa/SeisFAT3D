# ifndef STREAMER_HPP
# define STREAMER_HPP

# include "../geometry.hpp"

class Streamer : public Geometry
{
private:

    void set_relational();

public:

    void set_geometry(std::string file);

};

# endif
