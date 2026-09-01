# ifndef MIGRATION_CUH
# define MIGRATION_CUH

# include "../modeling/eikonal_iso.cuh"
# include "../modeling/eikonal_ani.cuh"

# define EPS 1e-6f

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
    int src_csIdx, rec_csIdy, curr;

    float ds, dr, dt, da, dCMP;
    float minCMPy, minCMPx;
    float maxCMPy, maxCMPx;
    float fmax, max_angle;
    float max_offset;
    float CMPx, CMPy;
    float aperture;

    bool anisotropy;

    float * h_xsrc, * h_ysrc;
    float * d_xsrc, * d_ysrc;

    float * h_xrec[2], * h_yrec[2];
    float * d_xrec[2], * d_yrec[2];

    float * h_data[2];
    float * d_data[2];

    float * h_Ts, * h_Tr[2];
    float * d_Ts, * d_Tr[2];

    size_t volBytes;
    cudaStream_t stream_cpy, stream_krn;
    cudaEvent_t cpy_done[2], krn_done[2], src_cpy_done;

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
    std::string zpos, xpos, ypos;

    Modeling * modeling = nullptr;

    void set_interpolation();
    void set_CMP_gathers();
    void set_anisotropy();
    void set_slowness();
    void set_wavelet();

    void perform_cubic(float * input, float * output);

    void show_information();    
    void prepare_convolution();
    void get_src_travel_times();
    void get_rec_travel_times();

    void set_src_line_components();

    virtual void set_migration() = 0;
    virtual void perform_forward() = 0;
    virtual void perform_adjoint() = 0;

public:
    
    std::string parameters;

    void set_parameters();

    virtual void kirchhoff_depth_migration() = 0;

    virtual void export_outputs() = 0;
};

__global__ void image_domain_adjoint_kernel(const float * __restrict__ d_S, const float * __restrict__ d_T_src, const float * __restrict__ d_T_rec, const float * __restrict__ data,
                                            float * __restrict__ model, const float * cs_xsrc, const float * cs_ysrc, const float * cs_xrec, const float * cs_yrec, const float aperture,
                                            const float max_offset, const int sm_z, const int sm_x, const int sm_y, const float cubic_dh, const float dh, const float dt, const int cubic_nxx,
                                            const int cubic_nyy, const int cubic_nzz, const int cubic_nb, const int nx, const int ny, const int nz, const int nt, const int nsy, const int nrx);

__global__ void image_domain_forward_kernel();

__global__ void angle_domain_adjoint_kernel();

__global__ void angle_domain_forward_kernel();

__device__ __forceinline__ float cubic1d(float p0, float p1, float p2, float p3, float u); 

__device__ __forceinline__ float tricubic(const float * __restrict__ smem_buffer, int base_z, int base_x, int base_y, 
                                          int smem_z, int smem_x, float uz, float ux, float uy); 

typedef struct 
{
    float T, dT_dz, dT_dx, dT_dy;
    float d2T_dz2, d2T_dx2, d2T_dy2;
    float d2T_dxdz, d2T_dydz, d2T_dxdy;
} TTD;

__device__ __forceinline__ TTD compute_TTDs(const float * __restrict__ smem_buffer, int base_z, int base_x, int base_y, 
                                            int smem_z, int smem_x, float uz, float ux, float uy, float cubic_dh); 

__device__ __forceinline__ void get_catmull_weights(float u, float w[4], float dw[4], float d2w[4]); 



# endif