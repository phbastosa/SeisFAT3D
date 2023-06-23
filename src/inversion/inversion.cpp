# include "inversion.hpp"

void Inversion::export_results()
{
    modeling->get_runtime();


    export_binary_float("gradient.bin", gradient, n_model);
}