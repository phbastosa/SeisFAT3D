# include "../migration/lfreq_migration.cuh"
# include "../migration/hfreq_migration.cuh"

int main(int argc, char **argv)
{
    std::vector<Migration *> migration = {
        new hfreq_Migration(),
        new lfreq_Migration()
    };

    auto file = std::string(argv[1]);
    auto type = std::stoi(catch_parameter("migration_type", file));

    migration[type]->file = file;

    migration[type]->set_parameters();

    for (int shot = 0; shot < migration[type]->total_shots; shot++)
    {
        migration[type]->shot_index;

        migration[type]->source_to_receiver_propagation();
        migration[type]->receiver_to_source_propagation();
    }

    migration[type]->image_compensation();

    migration[type]->export_seismic_volume();

    return 0;
}