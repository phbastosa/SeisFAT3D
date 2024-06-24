# ifndef MIGRATION_HPP
# define MIGRATION_HPP

# include <cuda_runtime.h>

# include "../modeling/modeling.hpp"

class Migration
{
protected:

    float * image = nullptr;

    std::string output_image_file;
    std::string output_image_folder;


    Modeling * modeling = nullptr;

    virtual void set_modeling_type() = 0;

public:

    std::string file;    

    int shot_index;
    int total_shots;
    int total_nodes;

    void set_parameters();
    void print_information();
    void import_seismic_data();

    void source_to_receiver_propagation();
    void receiver_to_source_propagation();

    void image_compensation();
    void export_seismic_volume();

};

# endif