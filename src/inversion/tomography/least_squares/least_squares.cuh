# ifndef LEAST_SQUARES_CUH
# define LEAST_SQUARES_CUH

# include "../tomography.hpp"

class Least_Squares : public Tomography
{
private:


protected:


public:

    void optimization();
    void set_parameters();
    void forward_modeling();
};

# endif