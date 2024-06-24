# ifndef MIGRATION_HPP
# define MIGRATION_HPP

# include <cuda_runtime.h>

# include "../modeling/modeling.hpp"

class Migration
{
private:

    void print_information();

protected:

    float * image = nullptr;

    std::string type_name;
    std::string type_message;
    std::string output_image_file;
    std::string output_image_folder;


    Modeling * modeling = nullptr;

    virtual void set_modeling_type() = 0;



    void source_to_receiver_propagation();
    void receiver_to_source_propagation();

public:

    std::string file;    

    void set_parameters();
    void import_gathers();
    void cross_correlation();
    void image_compensation();
    void export_seismic_volume();
};

# endif