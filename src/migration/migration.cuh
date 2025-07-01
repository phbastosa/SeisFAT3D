# ifndef MIGRATION_CUH
# define MIGRATION_CUH

# include "../modeling/eikonal_iso.cuh"
# include "../modeling/eikonal_ani.cuh"

class Migration
{
private:

    int nt; 
    int nBlocks; 
    int nThreads;

    float dt; 
    float aperture;
    float max_offset;

    float * d_Tr = nullptr;

    float * f_image = nullptr;
    float * h_image = nullptr;
    float * d_image = nullptr;

    float * h_seismic = nullptr;
    float * d_seismic = nullptr;

    std::string input_data_folder;
    std::string input_data_prefix;

    std::string output_image_folder;
    std::string output_table_folder;

    void show_information();
    void read_seismic_data();
    void set_receiver_point();
    void get_receiver_eikonal();
    void run_cross_correlation();
    void export_receiver_eikonal();

protected:

    Modeling * modeling = nullptr;

    virtual void set_modeling_type() = 0;
    
public:
    
    std::string parameters;

    void set_parameters();
    void image_building();
    void export_outputs();
};

__global__ void cross_correlation(float * Ts, float * Tr, float * image, float * seismic, float aperture_x, float aperture_y, 
                                  float cmp_x, float cmp_y, int spread, int nx, int ny, int nz, float dx, float dy, float dz, 
                                  int snx, int sny, int snz, float sdx, float sdy, float sdz, float scale, int nt, float dt);

# endif