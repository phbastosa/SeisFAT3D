# ifndef INVERSION_HPP
# define INVERSION_HPP

# include <string>
# include <vector> 
# include <iostream>

class Inversion
{
protected:

    std::string name;        

public:
    
    void get_name();

    virtual void set_name() = 0;
};

# endif