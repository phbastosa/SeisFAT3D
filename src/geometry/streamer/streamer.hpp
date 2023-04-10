# ifndef STREAMER_HPP
# define STREAMER_HPP

# include "../geometry.hpp"

class Streamer : public Geometry
{


public:    
    
    void set_parameters(std::string file);
};

# endif