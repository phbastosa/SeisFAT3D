#!/bin/bash

# Input Output scripts --------------------------------------------------------------------------------

io="../src/io/io.cpp"

# Acquisition geometry scripts ------------------------------------------------------------------------

geometry="../src/geometry/geometry.cpp"

regular="../src/geometry/regular/regular.cpp"
circular="../src/geometry/circular/circular.cpp"

geometry_main="../src/main/geometry_main.cpp"

geometry_all="$geometry $regular $circular"

# Seismic modeling scripts ----------------------------------------------------------------------------

eikonal="../src/modeling/eikonal.cpp"

pod="../src/modeling/podvin_and_lecomte/podvin_and_lecomte.cu"
fim="../src/modeling/fast_iterative_method/block_FIM.cu"
fsm="../src/modeling/fast_sweeping_method/accurate_FSM.cu"
ifim="../src/modeling/improved_FIM/improved_FIM.cu"

modeling_main="../src/main/modeling_main.cpp"

modeling_all="$eikonal $pod $fim $fsm $ifim"

acoustic_main="../src/main/acoustic_main.cpp"
acoustic_class="../src/seismogram/acoustic.cu"

# Seismic inversion scripts ---------------------------------------------------------------------------

tomography="../src/inversion/tomography.cpp"

least_squares="../src/inversion/least_squares/least_squares.cu"
adjoint_state="../src/inversion/adjoint_state/adjoint_state.cu"

inversion_main="../src/main/inversion_main.cpp"

inversion_all="$tomography $least_squares $adjoint_state"

# Seismic migration scripts ---------------------------------------------------------------------------

migration="../src/migration/migration.cpp"

kirchhoff="../src/migration/kirchhoff/kirchhoff.cu"

migration_main="../src/main/migration_main.cpp"

migration_all="$kirchhoff $migration"

# Compiler flags --------------------------------------------------------------------------------------

flags="--std=c++11 -lm -O3 -w -g"

# Main dialogue ---------------------------------------------------------------------------------------

USER_MESSAGE="
Usage:\n
    $ $0 -compile             # Create executables 
    $ $0 -modeling            # Perform eikonal solver          
    $ $0 -inversion           # Perform first arrival tomography
    $ $0 -migration           # Perform kirchhoff depth migration   
    $ $0 -seismogram          # Perform acoustic wave propagation

Tests:\n
    $ $0 -test_modeling       # Perform a small modeling experiment          
    $ $0 -test_inversion      # Perform a small inversion experiment
    $ $0 -test_migration      # Perform a small migration experiment          

Tools:\n
    $ $0 -configuration       # Check initial configuration plot
    $ $0 -manual_picking      # Perform manual picking algorithm using .segy files 

Projects:\n
    $ $0 -project_EAGE_2024   # Perform and show experiments of expanded abstract
    $ $0 -project_IMAGE_2024  # Perform and show experiments of expanded abstract 
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

    echo -e "../bin/\033[31mgeometry.exe\033[m" 
    nvcc $io $geometry_all $geometry_main $flags -o ../bin/geometry.exe

    echo -e "../bin/\033[31mmodeling.exe\033[m" 
    nvcc $io $geometry_all $modeling_all $modeling_main $flags -o ../bin/modeling.exe

    echo -e "../bin/\033[31minversion.exe\033[m" 
    nvcc $io $geometry_all $modeling_all $inversion_all $inversion_main $flags -o ../bin/inversion.exe

    echo -e "../bin/\033[31mmigration.exe\033[m"
    nvcc $io $geometry_all $modeling_all $migration_all $migration_main $flags -o ../bin/migration.exe

    echo -e "../bin/\033[31macoustic.exe\033[m"
    nvcc $io $geometry_all $acoustic_class $acoustic_main $flags -o ../bin/acoustic.exe

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

-seismogram) 
    
    ./../bin/acoustic.exe parameters.txt
    python3 ../tools/picking/add_trace_header.py parameters.txt
	
    exit 0
;;

-configuration)

    ./../bin/geometry.exe parameters.txt
    python3 ../tools/visualization/check_configuration.py parameters.txt
    
    exit 0
;;

-manual_picking)

    python3 ../tools/picking/manual_picking.py parameters.txt

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

	exit 0
;;

-test_inversion) 

    python3 ../tests/inversion/generate_models.py

    # ./../bin/modeling.exe ../tests/inversion/parFiles/parameters_obsData.txt

    # ./../bin/inversion.exe ../tests/inversion/parFiles/parameters_leastSquares.txt
    # ./../bin/inversion.exe ../tests/inversion/parFiles/parameters_adjointState.txt

    # ./../bin/modeling.exe ../tests/inversion/parFiles/parameters_leastSquares_finalModeling.txt
    # ./../bin/modeling.exe ../tests/inversion/parFiles/parameters_adjointState_finalModeling.txt

    python3 ../tests/inversion/generate_figures.py
	
    exit 0
;;

-test_migration)

    echo "testing a small migration experiment"

	exit 0
;;

-project_EAGE_2024)

    python3 ../projects/EAGE_2024/generate_models.py

    ./../bin/modeling.exe ../projects/EAGE_2024/parFiles/modFIM_parameters.txt
    ./../bin/modeling.exe ../projects/EAGE_2024/parFiles/modFSM_parameters.txt
    ./../bin/modeling.exe ../projects/EAGE_2024/parFiles/modIFIM_parameters.txt

    python3 ../projects/EAGE_2024/generate_figures.py

	exit 0
;;

-project_IMAGE_2024)



	exit 0
;;

* ) 

	echo -e "\033[31mERRO: Option $1 unknown!\033[m"
	echo -e "\033[31mType $0 -h for help \033[m"
	
    exit 3
;;

esac
