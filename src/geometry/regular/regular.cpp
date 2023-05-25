# include "regular.hpp"

void Regular::set_geometry()
{
    Geometry::set_geometry();
    
    if (import_geometry) 
    {
        import_coordinates();
    }
    else
    {
        std::vector<std::string> names = {"shots", "nodes"};
    
        for (auto name : names)
        {
            splitted = split(catch_parameter(name + "_nlines", file), ',');
            for (auto key : splitted) nlines.push_back(std::stoi(key));

            splitted = split(catch_parameter(name + "_SW", file), ',');
            for (auto key : splitted) SW.push_back(std::stof(key));

            splitted = split(catch_parameter(name + "_NW", file), ',');
            for (auto key : splitted) NW.push_back(std::stof(key));

            splitted = split(catch_parameter(name + "_SE", file), ',');
            for (auto key : splitted) SE.push_back(std::stof(key));

            if (name == std::string("shots")) set_regular(shots);
            if (name == std::string("nodes")) set_regular(nodes);
        }

        if (reciprocity) 
            set_reciprocity();

        export_coordinates();
    }   
}

