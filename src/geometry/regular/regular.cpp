# include "regular.hpp"

void Regular::set_parameters(std::string file)
{
    std::vector<std::string> splitted;

    name = "Regular geometry is active!";

    reciprocity = str2bool(catch_parameter("reciprocity", file));
    
    import_geometry = str2bool(catch_parameter("import_geometry", file));

    geometry_folder = catch_parameter("geometry_folder", file);

    splitted = split(catch_parameter("shots_SW", file), ',');
    for (auto key : splitted) shots_SW.push_back(std::stof(key));

    splitted = split(catch_parameter("shots_NW", file), ',');
    for (auto key : splitted) shots_NW.push_back(std::stof(key));

    splitted = split(catch_parameter("shots_SE", file), ',');
    for (auto key : splitted) shots_SE.push_back(std::stof(key));

    splitted = split(catch_parameter("nodes_SE", file), ',');
    for (auto key : splitted) nodes_SE.push_back(std::stof(key));

    splitted = split(catch_parameter("nodes_SE", file), ',');
    for (auto key : splitted) nodes_SE.push_back(std::stof(key));

    splitted = split(catch_parameter("nodes_SE", file), ',');
    for (auto key : splitted) nodes_SE.push_back(std::stof(key));

    splitted = split(catch_parameter("shots_nlines", file), ',');
    for (auto key : splitted) shots_nlines.push_back(std::stoi(key));

    splitted = split(catch_parameter("nodes_nlines", file), ',');
    for (auto key : splitted) nodes_nlines.push_back(std::stoi(key));

    import_geometry = str2bool(catch_parameter("import_geometry", file));
    geometry_folder = catch_parameter("geometry_folder", file);
    

}

void Regular::build_shots()
{



}

void Regular::build_nodes()
{




}