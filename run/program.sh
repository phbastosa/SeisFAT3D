#!/bin/bash

geometry="../src/geometry/geometry.cpp"
regular="../src/geometry/regular/regular.cpp"
circular="../src/geometry/circular/circular.cpp"
streamer="../src/geometry/streamer/streamer.cpp"
geometry_main="../src/geometry_main.cpp"

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
	echo " "
	echo "You didn't provide any parameter!" 
	echo "Type $0 -help for more info"
	echo " "  
	
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

    echo -e "../bin/\033[31minversion.exe\033[m" 

    echo -e "../bin/\033[31mmigration.exe\033[m" 


	exit 0
;;

-geometry) 

	exit 0
;;

-modeling) 
    
	exit 0
;;

-inversion) 
    
	exit 0
;;

-migration) 
    
	exit 0
;;

* ) ## Message for bad parameter
	
	echo -e "\033[31mERRO: Option $1 unknown!\033[m"
	echo -e "\033[31mType $0 -h for help \033[m"
	
    exit 3
;;

esac
