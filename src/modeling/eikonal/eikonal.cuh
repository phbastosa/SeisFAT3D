# ifndef EIKONAL_CUH
# define EIKONAL_CUH

# include "../modeling.hpp"

class Eikonal : public Modeling
{
protected:

    float * T = nullptr;
    
    int nit, nxx, nyy, nzz, volsize;

    void get_travelTimes();
    void get_firstArrivals();

    virtual void expansion() = 0;
    virtual void reduction() = 0;
    virtual void parameters() = 0;
    virtual void components() = 0;

public:

    void set_components();
    void build_outputs();    

    void set_parameters(std::string file);
};

# endif
