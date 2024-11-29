#!/bin/bash

# Input Output scripts --------------------------------------------------------------------------------

ioFunctions="../src/ioFunctions/ioFunctions.cpp"

# Acquisition geometry scripts ------------------------------------------------------------------------

geometry="../src/geometry/geometry.cpp"

# Seismic modeling scripts ----------------------------------------------------------------------------

modeling="../src/modeling/modeling.cpp"

eikonal="../src/modeling/hfreq/eikonal.cpp"
elastic="../src/modeling/lfreq/elastic.cu"

eikonal_iso="../src/modeling/hfreq/eikonal_iso.cu"
elastic_iso="../src/modeling/lfreq/elastic_iso.cu"

modeling_main="../src/modeling_main.cpp"

modeling_all="$modeling $eikonal $elastic $eikonal_iso $elastic_iso"

# Seismic inversion scripts ---------------------------------------------------------------------------

tomography="../src/inversion/tomography.cpp"

least_squares="../src/inversion/least_squares.cpp"
adjoint_state="../src/inversion/adjoint_state.cu"

inversion_main="../src/inversion_main.cpp"

inversion_all="$tomography $least_squares $adjoint_state"

# Seismic migration scripts ---------------------------------------------------------------------------

kirchhoff="../src/migration/kirchhoff.cu"

migration="../src/migration/migration.cpp"

migration_main="../src/migration_main.cpp"

migration_all="$migration $kirchhoff"

# Compiler flags --------------------------------------------------------------------------------------

flags="-Xcompiler -fopenmp --std=c++11 -lm -lfftw3 -O3"

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
    nvcc $ioFunctions $geometry $modeling_all $modeling_main $flags -o ../bin/modeling.exe

    echo -e "../bin/\033[31minversion.exe\033[m" 
    nvcc $ioFunctions $geometry $modeling_all $inversion_all $inversion_main $flags -o ../bin/inversion.exe

    echo -e "../bin/\033[31mmigration.exe\033[m"
    nvcc $ioFunctions $geometry $modeling_all $migration_all $migration_main $flags -o ../bin/migration.exe

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

    python3 -B ../tests/modeling/generate_models.py
    python3 -B ../tests/modeling/generate_geometry.py

    ./../bin/modeling.exe ../tests/modeling/parameters_eikonal.txt
    ./../bin/modeling.exe ../tests/modeling/parameters_elastic.txt

    python3 -B ../tests/modeling/generate_figures.py

	exit 0
;;

-test_inversion) 

    # python3 -B ../tests/inversion/generate_models.py
    # python3 -B ../tests/inversion/generate_geometry.py

    # ./../bin/modeling.exe ../tests/inversion/parameters_obsData.txt

    # ./../bin/inversion.exe ../tests/inversion/parameters_least_squares.txt
    ./../bin/inversion.exe ../tests/inversion/parameters_adjoint_state.txt

    python3 -B ../tests/inversion/generate_figures.py

    exit 0
;;

-test_migration)

    python3 -B ../tests/migration/generate_models.py
    python3 -B ../tests/migration/generate_geometry.py

    ./../bin/modeling.exe ../tests/migration/parameters.txt

    python3 -B ../tests/migration/data_preconditioning.py

    ./../bin/migration.exe ../tests/migration/parameters.txt

    python3 -B ../tests/migration/generate_figures.py

	exit 0
;;

* ) 

	echo -e "\033[31mERRO: Option $1 unknown!\033[m"
	echo -e "\033[31mType $0 -h for help \033[m"
	
    exit 3
;;

esac
