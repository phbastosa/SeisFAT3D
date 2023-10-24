#!/bin/bash

# Input Output scripts --------------------------------------------------------------------------------

io="../src/io/io.cpp"

# Acquisition geometry scripts ------------------------------------------------------------------------

geometry="../src/geometry/geometry.cpp"
regular="../src/geometry/regular/regular.cpp"
circular="../src/geometry/circular/circular.cpp"

geometry_all="$geometry $regular $circular"

# Seismic modeling scripts ----------------------------------------------------------------------------

modeling="../src/modeling/modeling.cpp"

pod="../src/modeling/podvin_and_lecomte/podvin_and_lecomte.cu"
fim="../src/modeling/fast_iterative_method/block_FIM.cu"
fsm="../src/modeling/fast_sweeping_method/accurate_FSM.cu"

modeling_main="../src/main/modeling_main.cpp"

modeling_all="$modeling $pod $fim $fsm $modeling_main"

# Seismic inversion scripts ---------------------------------------------------------------------------

inversion="../src/inversion/inversion.cpp"

tomography="../src/inversion/tomography/tomography.cpp"
least_squares="../src/inversion/tomography/least_squares/least_squares.cu"
adjoint_state="../src/inversion/tomography/adjoint_state/adjoint_state.cu"

inversion_main="../src/main/inversion_main.cpp"

inversion_all="$inversion $tomography $least_squares $adjoint_state $inversion_main"

# Seismic migration scripts ---------------------------------------------------------------------------

migration="../src/migration/migration.cpp"
kirchhoff="../src/migration/kirchhoff/kirchhoff.cpp"

migration_main="../src/main/migration_main.cpp"

migration_all="$migration $kirchhoff $migration_main"

# Compiler flags --------------------------------------------------------------------------------------

flags="--std=c++11 -lm -O3"

# Main dialogue ---------------------------------------------------------------------------------------

USER_MESSAGE="
Usage:\n
    $ $0 -compile        # Create executables 
    $ $0 -modeling       # Perform eikonal solver          
    $ $0 -inversion      # Perform first arrival tomography
    $ $0 -migration      # Perform kirchhoff depth migration   

Tests:\n
    $ $0 -test_modeling       # Perform a small modeling experiment          
    $ $0 -test_inversion      # Perform a small inversion experiment
    $ $0 -test_migration      # Perform a small migration experiment          
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

    echo -e "Compiling the stand-alone executables!\n"

    echo -e "../bin/\033[31mmodeling.exe\033[m" 
    nvcc $io $geometry_all $modeling_all $flags -o ../bin/modeling.exe

    # echo -e "../bin/\033[31minversion.exe\033[m" 
    # nvcc $io $geometry_all $modeling_all $inversion_all $flags -o ../bin/inversion.exe

    # echo -e "../bin/\033[31mmigration.exe\033[m"
    # nvcc $io $geometry_all $modeling_all $migration_all $flags -o ../bin/migration.exe

	exit 0
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

    python3 ../tests/modeling/generate_models.py

    spacings=(100 50 25)
    methods=("pod" "fim" "fsm")

    for method in ${methods[@]}; do 
        for spacing in ${spacings[@]}; do 
            ./../bin/modeling.exe ../tests/modeling/parFiles/parameters_"$method"_"$spacing"m.txt; 
        done    
    done 

    python3 ../tests/modeling/generate_figures.py

    rm runTime.txt

	exit 0
;;

-test_inversion)

    echo "testing a small inversion experiment"

	exit 0
;;

-test_migration)

    echo "testing a small migration experiment"

	exit 0
;;

* ) 

	echo -e "\033[31mERRO: Option $1 unknown!\033[m"
	echo -e "\033[31mType $0 -h for help \033[m"
	
    exit 3
;;

esac
