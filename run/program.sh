#!/bin/bash

# Input Output functions ---------------------------------------------------------------------------------

io="../src/utils/input_output/io.cpp"

# Acquisition geometry functions ------------------------------------------------------------------------

geometry="../src/geometry/geometry.cpp"
regular="../src/geometry/regular/regular.cpp"
circular="../src/geometry/circular/circular.cpp"

geometry_all="$geometry $regular $circular"

# Seismic modeling functions -----------------------------------------------------------------------------

modeling="../src/modeling/modeling.cpp"

eikonal="../src/modeling/eikonal_equation/eikonal.cpp"
pod="../src/modeling/eikonal_equation/podvin_and_lecomte/podvin_and_lecomte.cu"
fsm="../src/modeling/eikonal_equation/fast_sweeping_method/fast_sweeping_method.cu"
fim="../src/modeling/eikonal_equation/fast_iterative_method/fast_iterative_method.cu"

eikonal_all="$eikonal $pod $fim $fsm"

wave="../src/modeling/wave_equation/wave.cu"
scalar="../src/modeling/wave_equation/scalar/scalar.cu"
acoustic="../src/modeling/wave_equation/acoustic/acoustic.cu"
elastic="../src/modeling/wave_equation/elastic/elastic.cu"

modeling_main="../src/main/modeling_main.cpp"

wave_all="$wave $scalar $acoustic $elastic"

# Seismic inversion functions ----------------------------------------------------------------------------

inversion="../src/inversion/inversion.cpp"

tomography="../src/inversion/tomography/tomography.cpp"
least_squares="../src/inversion/tomography/least_squares/least_squares.cu"
adjoint_state="../src/inversion/tomography/adjoint_state/adjoint_state.cu"

tomography_all="$tomography $least_squares $adjoint_state"

waveform="../src/inversion/waveform/waveform.cpp"
# scalar_iso_fwi="../src/inversion/waveform/scalar_isotropic_fwi/scalar_isotropic_fwi.cu"

waveform_all="$waveform $scalar_fwi"

inversion_main="../src/main/inversion_main.cpp"

# Seismic migration functions ----------------------------------------------------------------------------

# migration="../src/migration/migration.cpp"
# kirchhoff="../src/migration/kirchhoff/kirchhoff.cpp"
# migration_main="../src/main/migration_main.cpp"

# Path unification ---------------------------------------------------------------------------------------

modeling_all="$modeling $eikonal_all $wave_all"
inversion_all="$inversion $tomography_all $waveform_all"
migration_all=""

flags="-Xcompiler=-fopenmp --std=c++11 --relocatable-device-code=true -lm -O3"

# Main dialogue ------------------------------------------------------------------------------------------

USER_MESSAGE="
Usage:\n
    $ $0 -compile        # Create executables 
    $ $0 -modeling       # Perform eikonal solver          
    $ $0 -inversion      # Perform adjoint-state tomography
    $ $0 -migration      # Perform kirchhoff depth migration   
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

    # echo -e "../bin/\033[31mmodeling.exe\033[m" 
    # nvcc $io $geometry_all $modeling_all $modeling_main $flags -o ../bin/modeling.exe

    echo -e "../bin/\033[31minversion.exe\033[m" 
    nvcc $io $geometry_all $modeling_all $inversion_all $inversion_main $flags -o ../bin/inversion.exe

    # echo -e "../bin/\033[31mmigration.exe\033[m"
    # nvcc $io $geometry $regular $circular $modeling $eikonal $migration $kirchhoff $migration_main $flags -o ../bin/migration.exe

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

* ) 

	echo -e "\033[31mERRO: Option $1 unknown!\033[m"
	echo -e "\033[31mType $0 -h for help \033[m"
	
    exit 3
;;

esac
