# include "lfreq_migration.cuh"

void lfreq_Migration::set_modeling_type()
{
    modeling = new lfreq_Modeling();

    type_message = "[1] Prestack Acoustic Reverse Time migration";
}