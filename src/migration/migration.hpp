# ifndef MIGRATION_HPP
# define MIGRATION_HPP

# include <string>
# include <vector> 
# include <iostream>

class Migration
{
protected:

    std::string name;        
    virtual void set_name() = 0;

public:
    
    void get_name();
};

# endif