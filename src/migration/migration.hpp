# ifndef MIGRATION_HPP
# define MIGRATION_HPP

# include "../modeling/modeling.hpp"

class Migration
{
private:

    std::string input_data_folder;
    std::string input_data_prefix;

protected:

    int nt;
    int nr;
    int ns;

    float dt;

    float * Tr = nullptr;
    float * Ts = nullptr;
    float * Im = nullptr;

    float * dTx = nullptr;
    float * dTy = nullptr;
    float * dTz = nullptr;

    float * data = nullptr;

    float * image = nullptr;

    Modeling * modeling = nullptr;

    virtual void set_modeling_type() = 0;

public:
    
    std::string file;

    void set_parameters();
    void read_input_data();

    virtual void image_building() = 0;

};

# endif