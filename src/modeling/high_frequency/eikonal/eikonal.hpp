# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include "../../modeling.hpp"

class Eikonal : public Modeling
{
private:


protected:

    float t0;    

    float * V = nullptr;
    float * S = nullptr;
    float * T = nullptr;

    void get_travel_times();
    void get_first_arrivals();
    void set_specifications();

    virtual void set_model_boundaries() = 0;
    virtual void set_modeling_message() = 0;
    virtual void set_preconditioners() = 0;

public:

    void build_outputs();

    virtual void initial_setup() = 0;
    virtual void forward_solver() = 0;
    virtual void free_space() = 0;
};

# endif