# ifndef KIRCHHOFF_CUH
# define KIRCHHOFF_CUH

# include "../migration.hpp"

class Kirchhoff : public Migration
{
protected:

    std::string name;        

public:
    
    void set_name();
};

# endif