#!/bin/bash

admin="../src/admin/admin.cpp"
geometry="../src/geometry/geometry.cpp"

# Seismic modeling scripts ----------------------------------------------------------------------------

folder="../src/modeling"

modeling="$folder/modeling.cu"

eikonal_iso="$folder/eikonal_iso.cu"
eikonal_ani="$folder/eikonal_ani.cu"

modeling_main="../src/modeling_main.cpp"

modeling_all="$modeling $eikonal $eikonal_iso $eikonal_ani"

# Seismic inversion scripts ---------------------------------------------------------------------------

folder="../src/inversion"

tomography="$folder/tomography.cpp"

least_squares="$folder/least_squares.cpp"
adjoint_state="$folder/adjoint_state.cu"

inversion_main="../src/inversion_main.cpp"

inversion_all="$tomography $least_squares $adjoint_state"

# Seismic migration scripts ---------------------------------------------------------------------------

folder="../src/migration"

kirchhoff="$folder/kirchhoff.cu"

migration="$folder/migration.cpp"

migration_main="../src/migration_main.cpp"

migration_all="$migration $kirchhoff"

# Compiler flags --------------------------------------------------------------------------------------

flags="-Xcompiler -fopenmp --std=c++11 --use_fast_math --relocatable-device-code=true -lm -O3"

# Main dialogue ---------------------------------------------------------------------------------------

USER_MESSAGE="
-------------------------------------------------------------------------------
                                 \033[34mSeisFAT2D\033[0;0m
-------------------------------------------------------------------------------
\nUsage:\n
    $ $0 -compile              
    $ $0 -modeling                      
    $ $0 -inversion           
    $ $0 -migration

Tests:

    $ $0 -test_modeling                      
    $ $0 -test_inversion           
    $ $0 -test_migration
    
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

    folder=../tests/modeling
    parameters=$folder/parameters.txt

    python3 -B $folder/generate_models.py
    python3 -B $folder/generate_geometry.py

    ./../bin/modeling.exe $parameters

    sed -i "s|modeling_type = 0|modeling_type = 1|g" "$parameters"

    ./../bin/modeling.exe $parameters

    sed -i "s|modeling_type = 1|modeling_type = 0|g" "$parameters"

    python3 -B $folder/generate_figures.py

	exit 0
;;

-test_inversion) 

    folder=../tests/inversion

    python3 -B $folder/generate_models.py
    python3 -B $folder/generate_geometry.py

    ./../bin/modeling.exe $folder/parameters_obsData.txt

    ./../bin/inversion.exe $folder/parameters_least_squares.txt
    ./../bin/inversion.exe $folder/parameters_adjoint_state.txt

    python3 -B $folder/generate_figures.py

    exit 0
;;

-test_migration)

    folder=../tests/migration

    python3 -B $folder/generate_models.py
    python3 -B $folder/generate_geometry.py

    ./../bin/modeling.exe $folder/parameters_obsData.txt

    python3 -B $folder/data_preconditioning.py

    ./../bin/migration.exe $folder/parameters_kirchhoff.txt

    python3 -B $folder/generate_figures.py

	exit 0
;;

* ) 

	echo -e "\033[31mERRO: Option $1 unknown!\033[m"
	echo -e "\033[31mType $0 -h for help \033[m"
	
    exit 3
;;

esac
