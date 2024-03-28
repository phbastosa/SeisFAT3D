# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include "../modeling.hpp"

class Eikonal : public Modeling
{
private:

    void set_models();

protected:


    virtual void set_volumes() = 0;
    virtual void set_specifics() = 0;
    virtual void initialization() = 0;

public: 

    virtual void set_forward_solver() = 0;

};

# endif