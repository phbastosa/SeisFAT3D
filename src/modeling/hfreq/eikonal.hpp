# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include "../modeling.hpp"

class Eikonal : public Modeling
{
private:

    void set_boundaries();
    void set_specifications();

protected:

    void compute_seismogram();

    virtual void set_conditions() = 0;
    virtual void set_properties() = 0;

public:

    void initialization();

    virtual void forward_solver() = 0;

    void export_synthetic_data();
};

# endif