# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include "../modeling.hpp"

class Eikonal : public Modeling
{
private:

    void set_models();
    void set_outputs();

protected:

    virtual void set_volumes() = 0;
    virtual void set_specifics() = 0;
    virtual void initialization() = 0;

    void get_receiver_output();
    void get_wavefield_output();

public: 

    virtual void set_forward_solver() = 0;
    virtual void free_space() = 0;
};

# endif