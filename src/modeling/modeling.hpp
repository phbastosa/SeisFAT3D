# ifndef MODELING_HPP
# define MODELING_HPP

# include <string>
# include <vector> 
# include <iostream>

class Modeling
{
protected:

    std::string name;        

public:
    
    void get_name();

    virtual void set_name() = 0;
};

# endif