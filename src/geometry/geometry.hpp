# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include <string>
# include <vector> 
# include <iostream>

class Geometry
{
protected:

    std::string name;        

public:
    
    void get_name();

    virtual void set_name() = 0;
};

# endif