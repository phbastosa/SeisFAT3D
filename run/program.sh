#!/bin/bash

io="../src/utils/input_output/io.cpp"

geometry="../src/geometry/geometry.cpp"
regular="../src/geometry/regular/regular.cpp"
circular="../src/geometry/circular/circular.cpp"

modeling="../src/modeling/modeling.cpp"

PAL="../src/modeling/eikonal/PAL.cu"
FSM="../src/modeling/eikonal/FSM.cu"
FIM="../src/modeling/eikonal/FIM.cu"
eikonal="../src/modeling/eikonal/eikonal.cpp"

scalar="../src/modeling/scalar/scalar.cpp"
acoustic="../src/modeling/acoustic/acoustic.cpp"
elastic="../src/modeling/elastic/elastic.cpp"
modeling_main="../src/main/modeling_main.cpp"

inversion="../src/inversion/inversion.cpp"
waveform="../src/inversion/waveform/waveform.cpp"
tomography="../src/inversion/tomography/tomography.cpp"
inversion_main="../src/main/inversion_main.cpp"

migration="../src/migration/migration.cpp"
kirchhoff="../src/migration/kirchhoff/kirchhoff.cpp"
reverseTime="../src/migration/reverseTime/reverseTime.cpp"
migration_main="../src/main/migration_main.cpp"

flags="-std=c++11 -lm -O3"

USER_MESSAGE="
Usage:
    $ $0 -help           # 
    $ $0 -compile        # 
    $ $0 -modeling       #           
    $ $0 -inversion      # 
    $ $0 -migration      #    
"

[ -z "$1" ] && 
{
	echo -e "\nYou didn't provide any parameter!" 
	echo -e "Type $0 -help for more info\n"
    exit 1 
}

case "$1" in

-help) 

	echo -e "$USER_MESSAGE"
	exit 0
;;

-compile) 

    echo -e "Compiling the stand-alone executables!\n"

    echo -e "../bin/\033[31mmodeling.exe\033[m" 
    nvcc $io $geometry $regular $circular $modeling $PAL $FIM $FSM $eikonal $scalar $acoustic $elastic $modeling_main $flags -o ../bin/modeling.exe

    # echo -e "../bin/\033[31minversion.exe\033[m" 
    # nvcc $io $geometry $regular $circular $modeling $PAL $FIM $FSM $eikonal $inversion $waveform $tomography $inversion_main $flags -o ../bin/inversion.exe

    # echo -e "../bin/\033[31mmigration.exe\033[m"
    # nvcc $io $migration $kirchhoff $reverseTime $migration_main -lm -O3 -o ../bin/migration.exe

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
