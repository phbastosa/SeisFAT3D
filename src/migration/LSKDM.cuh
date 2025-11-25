# ifndef LSKDM_CUH
# define LSKDM_CUH

# include "migration.cuh"

class LSKDM : public Migration
{
private:

    int iteration;

    bool converged;

    std::vector<float> residuo;

    float residuals;

    float beta, beta_num, beta_den;
    float alpha, alpha_num, alpha_den;
    
    void update_model();
    void initialization();
    void compute_gradient();
    void compute_residuals();
    void compute_direction();
    void compute_stepLength();
    void show_iteration_info();

protected:

    float * h_gradient = nullptr;
    float * d_gradient = nullptr;

    float * gradient_old = nullptr;
    
    float * h_direction = nullptr;
    float * d_direction = nullptr;

    virtual void set_migration() = 0;
    virtual void perform_forward() = 0;
    virtual void perform_adjoint() = 0;
    
    virtual void perform_adjoint_gradient() = 0;
    virtual void perform_forward_direction() = 0;

public:

    void kirchhoff_depth_migration();

    void export_outputs();
};

# endif
