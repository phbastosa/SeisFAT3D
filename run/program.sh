#!/bin/bash

# Input Output scripts --------------------------------------------------------------------------------

admin="../src/admin/admin.cpp"

# Acquisition geometry scripts ------------------------------------------------------------------------

geometry="../src/geometry/geometry.cpp"

# Seismic modeling scripts ----------------------------------------------------------------------------

modeling="../src/modeling/modeling.cu"

eikonal_iso="../src/modeling/eikonal_iso.cu"
eikonal_ani="../src/modeling/eikonal_ani.cu"

modeling_main="../src/modeling_main.cpp"

modeling_all="$modeling $eikonal_iso $eikonal_ani"

# Seismic inversion scripts ---------------------------------------------------------------------------

inversion="../src/inversion/inversion.cpp"

tomography_iso="../src/inversion/tomography_iso.cpp"
tomography_vti="../src/inversion/tomography_vti.cpp"

inversion_main="../src/inversion_main.cpp"

inversion_all="$inversion $tomography_iso $tomography_vti"

# Seismic migration scripts ---------------------------------------------------------------------------

migration="../src/migration/migration.cu"

kirchhoff_iso="../src/migration/kirchhoff_iso.cu"
kirchhoff_ani="../src/migration/kirchhoff_ani.cu"

migration_main="../src/migration_main.cpp"

migration_all="$migration $kirchhoff_iso $kirchhoff_ani"

# Compiler flags --------------------------------------------------------------------------------------

flags="-Xcompiler -fopenmp --std=c++11 --relocatable-device-code=true -lm -O3"

# Main dialogue ---------------------------------------------------------------------------------------

USER_MESSAGE="
-------------------------------------------------------------------------------
                                 \033[34mSeisFAT3D\033[0;0m
-------------------------------------------------------------------------------
\nUsage:
        $ $0 -compile              
        $ $0 -modeling                      
        $ $0 -inversion           
        $ $0 -migration
-------------------------------------------------------------------------------
"

[ -z "$1" ] && 
{
	echo -e "\nYou didn't provide any parameter!" 
	echo -e "Type $0 -help for more info\n"
    exit 1 
}

case "$1" in

-h) 

	echo -e "$USER_MESSAGE"
	exit 0
;;

-compile) 

    echo -e "Compiling stand-alone executables!\n"

    echo -e "../bin/\033[31mmodeling.exe\033[m" 
    nvcc $admin $geometry $modeling_all $modeling_main $flags -o ../bin/modeling.exe

    # echo -e "../bin/\033[31minversion.exe\033[m" 
    # nvcc $admin $geometry $modeling_all $inversion_all $inversion_main $flags -o ../bin/inversion.exe

    # echo -e "../bin/\033[31mmigration.exe\033[m"
    # nvcc $admin $geometry $modeling_all $migration_all $migration_main $flags -o ../bin/migration.exe

	exit 0
;;

-clean)

    rm ../bin/*.exe
    rm ../inputs/data/*.bin
    rm ../inputs/geometry/*.txt
    rm ../inputs/models/*.bin
    rm ../outputs/convergence/*.txt
    rm ../outputs/migratedImages/*.bin
    rm ../outputs/recoveredModels/*.bin
    rm ../outputs/syntheticData/*.bin
    rm ../outputs/travelTimeTables/*.bin
;;

-modeling) 

    ./../bin/modeling.exe parameters.txt
	
    exit 0
;;

-inversion) 
    
    ./../bin/inversion.exe parameters.txt
	
    exit 0
;;

-migration) 
    
    ./../bin/migration.exe parameters.txt
	
    exit 0
;;


-test_modeling)

    prefix=../tests/modeling
    parameters=$prefix/parameters.txt

    python3 -B $prefix/generate_models.py $parameters
    python3 -B $prefix/generate_geometry.py $parameters

    ./../bin/modeling.exe $parameters

    python3 -B $prefix/generate_figures.py $parameters

	exit 0
;;

-test_inversion) 

    prefix=../tests/inversion
    parameters=$prefix/parameters.txt

    python3 -B $prefix/generate_models.py $parameters
    python3 -B $prefix/generate_geometry.py $parameters

    true_model="model_file = ../inputs/models/inversion_test_true_vp.bin"
    init_model="model_file = ../inputs/models/inversion_test_init_vp.bin"

    ./../bin/modeling.exe $parameters

    sed -i "s|$true_model|$init_model|g" "$parameters"    
    
    ./../bin/inversion.exe $parameters

    sed -i "s|$init_model|$true_model|g" "$parameters"

    python3 -B $prefix/generate_figures.py $parameters

    exit 0
;;

-test_migration)

    prefix=../tests/migration
    parameters=$prefix/parameters.txt

    python3 -B $prefix/generate_models.py $parameters
    python3 -B $prefix/generate_geometry.py $parameters
    python3 -B $prefix/generate_input_data.py $parameters

    ./../bin/migration.exe $parameters

    python3 -B $prefix/generate_figures.py $parameters

	exit 0
;;

* ) 

	echo -e "\033[31mERRO: Option $1 unknown!\033[m"
	echo -e "\033[31mType $0 -h for help \033[m"
	
    exit 3
;;

esac