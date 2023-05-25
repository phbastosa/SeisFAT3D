# ifndef EIKONAL_FSM_CUH
# define EIKONAL_FSM_CUH

# include "eikonal.hpp"

class Eikonal_fsm : public Eikonal
{
    int i, j, k, padb;
    int sgnvx, sgnvy, sgnvz;
    int sgntx, sgnty, sgntz;

    float dxi, dyi, dzi;
    float dx2i, dy2i, dz2i, dsum;
    float dz2dx2, dz2dy2, dx2dy2;

    void expansion();
    void reduction();
    void parameters();
    void components();

    void init_sweep();
    void full_sweep();    
    void inner_sweep();

public:

    void free_space();
    void initial_setup();
    void forward_solver();
};

# endif