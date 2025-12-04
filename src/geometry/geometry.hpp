# ifndef GEOMETRY_HPP
# define GEOMETRY_HPP

# include "../admin/admin.hpp"

class Geometry
{
public:

    int nsrc;
    int nrec;

    float * xsrc = nullptr;
    float * ysrc = nullptr;
    float * zsrc = nullptr;

    float * xrec = nullptr;
    float * yrec = nullptr;
    float * zrec = nullptr;

    std::string parameters;

    void set_parameters();     
};

# endif