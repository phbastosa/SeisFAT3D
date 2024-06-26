# ifndef MODELING_HPP
# define MODELING_HPP

# include "../geometry/geometry.hpp"

class Modeling
{
private:

    std::chrono::system_clock::time_point ti, tf;

    void set_generals();
    void set_geometry();

    void check_geometry_overflow();
    
protected:

    std::string type_name;
    std::string type_message;

    bool export_receiver_output;
    int receiver_output_samples;

    std::string receiver_output_file;
    std::string receiver_output_folder;

    int sidx, sidy, sidz;

    void set_vp_model();
    void set_boundary();

    virtual void set_volumes() = 0;
    virtual void set_outputs() = 0;
    virtual void set_specifics() = 0;

    virtual void initialization() = 0;

    virtual void get_receiver_output() = 0;

public:

    std::string file;

    int shot_index;
    int time_index;
    int source_index;    

    int total_shots;
    int total_nodes;

    int blocksPerGrid;
    int threadsPerBlock;

    float dx, dy, dz;
    int nx, ny, nz, nPoints;
    int nxx, nyy, nzz, volsize;

    int nbxl, nbxr, nbyl; 
    int nbyr, nbzu, nbzd;

    Geometry * geometry;

    float * S = nullptr;
    float * V = nullptr;

    float * T = nullptr;
    float * P = nullptr;

    float * model = nullptr;

    float * receiver_output = nullptr;

    void expand_boundary(float * input, float * output);
    void reduce_boundary(float * input, float * output);

    void set_runtime();
    void get_runtime();

    void print_information();

    void set_parameters(); 
    
    void set_initial_conditions();
    
    void export_outputs();

    virtual void forward_propagation() = 0;
    virtual void free_space() = 0;    
};

# endif