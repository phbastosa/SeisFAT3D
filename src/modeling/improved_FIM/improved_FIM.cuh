# ifndef HD_FIM_HPP
# define HD_FIM_HPP

# include "../eikonal.hpp"

# include "travel_time_3d_3_meat.hpp"

# define BLOCK_LENGTH 4

using namespace jarvis;

class Improved_FIM : public Eikonal
{
private:

    travel_time_3d_module travel_time_3d_2rd_diag;

    void set_specific_boundary();
    void set_eikonal_volumes();
    void initialization();

public:

    void forward_solver();
    void free_space();
};

# endif