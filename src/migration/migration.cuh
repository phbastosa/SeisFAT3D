# ifndef MIGRATION_CUH
# define MIGRATION_CUH

# include "../modeling/eikonal_iso.cuh"
# include "../modeling/eikonal_ani.cuh"

# define EPS 1e-6f

# define NTHREADS 256

class Migration
{
protected:

    int cmpId, nCMP, nCMPx, nCMPy;
    int nt, nang, nw, nfft, max_it;
    int nBlocks, d_samples, m_samples; 

    float max_angle, da;
    float ds, dr, dt, fmax;
    float minCMPy, minCMPx;

    bool anisotropy;
    
    float * seismic = nullptr; 
    float * wavelet = nullptr;

    double * time_trace = nullptr;
    double * time_wavelet = nullptr;

    fftw_complex * freq_trace = nullptr;
    fftw_complex * freq_wavelet = nullptr;

    fftw_plan trace_forward_plan;
    fftw_plan trace_inverse_plan;
    fftw_plan wavelet_forward_plan;

    float * m1 = nullptr;
    float * d1 = nullptr;

    float * m2 = nullptr;
    float * d2 = nullptr;

    float * h_Ts = nullptr;
    float * h_Tr = nullptr;

    float * d_Ts = nullptr;
    float * d_Tr = nullptr;

    float * d_data = nullptr;
    float * h_data = nullptr;    

    float * h_model = nullptr;
    float * d_model = nullptr;

    std::string domain;
    std::string migType;
    std::string output_path;

    Modeling * modeling = nullptr;
    
    std::string current, keyword;
    std::string xpos, ypos, zpos, total;
    std::string current_operation;

    std::string input_data_folder;
    std::string input_data_prefix;
    
    std::string tables_folder;
    std::string images_folder;
    std::string gathers_folder;
    std::string residuo_folder;

    void set_wavelet();
    void set_gathers();

    void set_src_domain();
    void set_current_src();

    void set_rec_domain();
    void set_current_rec();

    void show_information();    
    void set_src_travel_times();
    void set_rec_travel_times();
    void prepare_convolution();

    void adjoint_convolution();
    void forward_convolution();

    virtual void set_migration() = 0;
    virtual void perform_forward() = 0;
    virtual void perform_adjoint() = 0;

public:
    
    std::string parameters;

    void set_parameters();

    void dot_product_test();

    virtual void kirchhoff_depth_migration() = 0;

    virtual void export_outputs() = 0;
};

__global__ void image_domain_forward_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dx, float dz, float dt, int nxx, int nzz, int nt, int nb);
__global__ void image_domain_adjoint_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dx, float dz, float dt, int nxx, int nzz, int nt, int nb);

__global__ void angle_domain_forward_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dx, float dz, float dt, float da, int nxx, int nzz, int nt, int na, int nb, int cmpId);
__global__ void angle_domain_adjoint_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dx, float dz, float dt, float da, int nxx, int nzz, int nt, int na, int nb, int cmpId);

# endif