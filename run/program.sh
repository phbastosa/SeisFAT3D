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
tomography_ort="../src/inversion/tomography_ort.cpp"

inversion_main="../src/inversion_main.cpp"

inversion_all="$inversion $tomography_iso"

# Seismic migration scripts ---------------------------------------------------------------------------

migration="../src/migration/migration.cu"

KDM="../src/migration/KDM.cu"
IDKDM="../src/migration/IDKDM.cu"
ADKDM="../src/migration/ADKDM.cu"

LSKDM="../src/migration/LSKDM.cu"
IDLSKDM="../src/migration/IDLSKDM.cu"
ADLSKDM="../src/migration/ADLSKDM.cu"

migration_main="../src/migration_main.cpp"

migration_all="$migration $KDM $LSKDM $IDKDM $IDLSKDM $ADKDM $ADLSKDM"

# Compiler flags --------------------------------------------------------------------------------------

flags="-Xcompiler -fopenmp -lfftw3 --std=c++11 --relocatable-device-code=true -lm -O3"

# Main dialogue ---------------------------------------------------------------------------------------

USER_MESSAGE="
-------------------------------------------------------------------------------
 \033[34mSeisFAT3D\033[0;0m --------------------------------------------------------------------
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

    echo -e "../bin/\033[31mmigration.exe\033[m"
    nvcc $admin $geometry $modeling_all $migration_all $migration_main $flags -o ../bin/migration.exe

	exit 0
;;

-clean)

    rm *.png
    rm ../bin/*.exe
    rm ../inputs/data/*.bin
    rm ../inputs/geometry/*.txt
    rm ../inputs/models/*.bin
    rm ../outputs/residuo/*.txt
    rm ../outputs/seismic/*.bin
    rm ../outputs/models/*.bin
    rm ../outputs/data/*.bin
    rm ../outputs/tables/*.bin
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