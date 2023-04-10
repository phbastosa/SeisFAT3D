# ifndef STREAMER_HPP
# define STREAMER_HPP

# include "../geometry.hpp"

class Streamer : public Geometry
{


public:    
    
    void set_geometry(std::string file, std::string name);
};

# endif