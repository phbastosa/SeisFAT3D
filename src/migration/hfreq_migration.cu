# include "hfreq_migration.cuh"

void hfreq_Migration::set_modeling_type()
{
    modeling = new hfreq_Modeling();

    type_message = "[0] Prestack Kirchhoff depth migration";
}