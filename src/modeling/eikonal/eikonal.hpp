# ifndef EIKONAL_HPP
# define EIKONAL_HPP

# include "../modeling.hpp"

class Eikonal : public Modeling
{
private:

    std::string vp_model_file;

    float * vp = nullptr;




public:

    float * T = nullptr;




    void set_parameters(std::string file);
};

# endif
