#!/bin/bash

geometry="../src/geometry/geometry.cpp"
regular="../src/geometry/regular/regular.cpp"
circular="../src/geometry/circular/circular.cpp"
streamer="../src/geometry/streamer/streamer.cpp"
geometry_main="../src/geometry_main.cpp"

modeling="../src/modeling/modeling.cpp"
eikonal="../src/modeling/eikonal/eikonal.cpp"
scalar="../src/modeling/scalar/scalar.cpp"
acoustic="../src/modeling/acoustic/acoustic.cpp"
elastic="../src/modeling/elastic/elastic.cpp"
modeling_main="../src/modeling_main.cpp"

inversion="../src/inversion/inversion.cpp"
waveform="../src/inversion/waveform/waveform.cpp"
tomography="../src/inversion/tomography/tomography.cpp"
inversion_main="../src/inversion_main.cpp"

migration="../src/migration/migration.cpp"
kirchhoff="../src/migration/kirchhoff/kirchhoff.cpp"
reverseTime="../src/migration/reverseTime/reverseTime.cpp"
migration_main="../src/migration_main.cpp"

USER_MESSAGE="
Usage:
    $ $0 -help           # 
    $ $0 -compile        # 
    $ $0 -geometry       #  
    $ $0 -modeling       #           
    $ $0 -inversion      # 
    $ $0 -migration      #    
"
# Check if user provide some parameter
[ -z "$1" ] && {
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

    echo -e "../bin/\033[31mgeometry.exe\033[m" 
    nvcc $geometry $regular $circular $streamer $geometry_main -lm -O3 -o ../bin/geometry.exe

    echo -e "../bin/\033[31mmodeling.exe\033[m" 
    nvcc $modeling $eikonal $scalar $acoustic $elastic $modeling_main -lm -O3 -o ../bin/modeling.exe

    echo -e "../bin/\033[31minversion.exe\033[m" 
    nvcc $inversion $waveform $tomography $inversion_main -lm -O3 -o ../bin/inversion.exe

    echo -e "../bin/\033[31mmigration.exe\033[m"
    nvcc $migration $kirchhoff $reverseTime $migration_main -lm -O3 -o ../bin/migration.exe

	exit 0
;;

-geometry) 

    ./../bin/geometry.exe

	exit 0
;;

-modeling) 

    ./../bin/modeling.exe

	exit 0
;;

-inversion) 
    
    ./../bin/inversion.exe

	exit 0
;;

-migration) 
    
    ./../bin/migration.exe

	exit 0
;;

* ) ## Message for bad parameter
	
	echo -e "\033[31mERRO: Option $1 unknown!\033[m"
	echo -e "\033[31mType $0 -h for help \033[m"
	
    exit 3
;;

esac
