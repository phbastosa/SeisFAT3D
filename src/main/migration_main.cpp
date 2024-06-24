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

    migration[type]->import_gathers();

    migration[type]->cross_correlation();

    migration[type]->image_compensation();

    migration[type]->export_seismic_volume();

    return 0;
}