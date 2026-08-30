# ifndef MIGRATION_CUH
# define MIGRATION_CUH

# include "../modeling/eikonal_iso.cuh"
# include "../modeling/eikonal_ani.cuh"

# define EPS 1e-6f

# define NTHREADS 256

class Migration
{
protected:

    int old_nx, old_ny, old_nz, old_nPoints;
    int new_nx, new_ny, new_nz, new_nPoints;

    float old_dh, new_dh;

    int nfreq, nsx, nsy, nrx, nry;
    int cmpId, nCMP, nCMPx, nCMPy;
    int nt, nw, nfft, max_it, nang;
    int nBlocks, d_samples, m_samples; 
    int src_csIdx, rec_csIdy;

    float ds, dr, dt, da, dCMP;
    float minCMPy, minCMPx;
    float maxCMPy, maxCMPx;
    float fmax, max_angle;
    float max_offset;
    float CMPx, CMPy;
    float aperture;

    bool anisotropy;

    float * h_data[2];
    float * d_data[2];

    float * h_Ts, * h_Tr[2];
    float * d_Ts, * d_Tr[2];

    size_t volBytes;
    cudaStream_t stream_cpy, stream_krn;
    cudaEvent_t cpy_done[2], krn_done[2];

    float * h_model = nullptr;
    float * d_model = nullptr;
    
    float * seismic = nullptr; 
    float * wavelet = nullptr;

    double * time_trace = nullptr;
    double * time_wavelet = nullptr;

    fftw_complex * freq_trace = nullptr;
    fftw_complex * freq_wavelet = nullptr;

    fftw_plan trace_forward_plan;
    fftw_plan trace_inverse_plan;
    fftw_plan wavelet_forward_plan;

    std::string domain;
    std::string modType;
    std::string migType;
    std::string output_path;
    
    std::string input_data_folder;
    std::string input_data_prefix;
    
    std::string tables_folder;
    std::string seismic_folder;
    std::string residuo_folder;

    Modeling * modeling = nullptr;

    void set_interpolation();
    void set_CMP_gathers();
    void set_anisotropy();
    void set_slowness();
    void set_wavelet();

    void perform_cubic(float * input, float * output);

    void show_information();    
    void prepare_convolution();
    void set_src_travel_times();
    void set_rec_travel_times();

    void set_src_line_seismic();
    void set_src_line_travel_times();

//    void adjoint_convolution();
//    void forward_convolution();

    virtual void set_migration() = 0;
    virtual void perform_forward() = 0;
    virtual void perform_adjoint() = 0;

public:
    
    std::string parameters;

    void set_parameters();

    //void dot_product_test();

    virtual void kirchhoff_depth_migration() = 0;

    virtual void export_outputs() = 0;
};

__global__ void image_domain_adjoint_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dt, int nt, 
                                            float old_dx, float old_dy, float old_dz, float new_dx, float new_dy, float new_dz, 
                                            int old_nx, int old_ny, int old_nz, int new_nxx, int new_nyy, int new_nzz, int nb,
                                            float aperture, float cmpx, float cmpy);

__global__ void image_domain_forward_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dt, int nt, 
                                            float old_dx, float old_dy, float old_dz, float new_dx, float new_dy, float new_dz, 
                                            int old_nx, int old_ny, int old_nz, int new_nxx, int new_nyy, int new_nzz, int nb,
                                            float aperture, float cmpx, float cmpy);

__global__ void angle_domain_adjoint_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dx, float dz, float dt, float da, int nxx, int nzz, int nt, int na, int nb, int cmpId);
__global__ void angle_domain_forward_kernel(float * S, float * Ts, float * Tr, float * data, float * model, float dx, float dz, float dt, float da, int nxx, int nzz, int nt, int na, int nb, int cmpId);

__device__ float d_cubic1d(float P[4], float dx);
__device__ float d_cubic2d(float P[4][4], float dx, float dy);
__device__ float d_cubic3d(float P[4][4][4], float dx, float dy, float dz);


# endif