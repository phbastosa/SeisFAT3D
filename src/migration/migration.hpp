# ifndef MIGRATION_HPP
# define MIGRATION_HPP

# include "../modeling/hfreq/eikonal_iso.cuh"

class Migration
{
private:

    std::string input_data_folder;
    std::string input_data_prefix;

    void initialization();
    void show_information();
    void get_receiver_traveltimes();
    void export_receiver_traveltimes();

protected:

    float * Tr = nullptr;
    float * Ts = nullptr;

    float * image = nullptr;
    float * seismic = nullptr;

    Modeling * modeling = nullptr;

    virtual void set_specifications() = 0;
    virtual void run_cross_correlation() = 0;

public:
    
    std::string parameters;

    void set_parameters();
    
    void read_seismic_data();

    void image_building();
    void export_outputs();
};

# endif