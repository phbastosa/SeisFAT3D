# ifndef KDM_CUH
# define KDM_CUH

# include "migration.cuh"

class KDM : public Migration
{
protected:
    
    virtual void set_migration() = 0;
    virtual void perform_forward() = 0;
    virtual void perform_adjoint() = 0;
    
public:

    void kirchhoff_depth_migration();

    void export_outputs();
};

# endif
