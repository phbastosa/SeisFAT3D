# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include "../modeling.hpp"

class Eikonal : public Modeling
{
protected:

    void set_model_boundaries();
    void set_slowness_model();
    void set_outputs();

    void get_travel_times();
    void get_first_arrivals();

public:

    float t0;    

    float * S = nullptr;
    float * T = nullptr;

    void build_outputs();

    virtual void set_parameters() = 0; 

    virtual void info_message() = 0;
    virtual void initial_setup() = 0;
    virtual void forward_solver() = 0;

    virtual void free_space() = 0;
};

# endif