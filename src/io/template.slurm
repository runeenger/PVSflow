
## Set up jo environment:
set -o errexit # Exit the script on any error
set -o nounset # Treat any unset variables as an error

#module --quiet purge # Reset the modules to the system default

module load matplotlib/3.0.0-intel-2018b-Python-3.6.6

## Set-up FEniCS 2019
source /cluster/shared/fenics/conf/fenics-2019.1.0.saga.intel.conf

source /cluster/projects/nn9279k/johannr/gmsh.conf
source /cluster/projects/nn9279k/johannr/ulfy.conf

# Work directory
workdir=$USERWORK/
mkdir -pv $workdir


# Copy to work directory

PROJECT_FOLDER=sleep
#yes|cp -ruv /cluster/home/avallet/${PROJECT_FOLDER} ${workdir}/ || echo "skip copying due to other thread"
yes|cp -ruv /cluster/projects/nn9279k/Alexandra/${PROJECT_FOLDER} ${workdir}/ || echo "skip copying due to other thread"



cd ${workdir}/${PROJECT_FOLDER}

# Add the project to the python path 
export PYTHONPATH=`pwd`:"$PYTHONPATH"

## Do work:
SCRIPT_FOLDER=sleep/fbb_DD/

cd ${SCRIPT_FOLDER}

echo 'launch the script'
