#!/bin/bash
#SBATCH --account=p31861  ## YOUR ACCOUNT pXXXX or bXXXX
#SBATCH --partition=short ### PARTITION (buyin, short, normal, etc)
#SBATCH --nodes=1 ## how many computers do you need
#SBATCH --ntasks-per-node=1 ## how many cpus or processors do you need on each computer
#SBATCH --time=03:59:00 ## how long does this need to run (remember different partitions have restrictions on this param)
#SBATCH --mem-per-cpu=3G ## how much RAM do you need per CPU (this effects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name=BMC  ## When you run squeue -u NETID this is how you can identify the job

# load modules
module purge
module load singularity
module load mesa

cd /projects/p31861/Users/LiuJiaqi/workdir/RobotGcode
rm -r build

mkdir -p /projects/p31861/Users/LiuJiaqi/workdir/RobotGcode/build

pwd

# configure chrono-concrete
singularity exec --pwd /RobotGcode/build -B /projects/p31861/Users/LiuJiaqi/workdir/RobotGcode:/RobotGcode -B /projects/p31861/Users/AkpanWisdom/NEWCHRONO/chrono-concrete:/chrono-concrete /projects/p31861/Users/LaleErol/project-chrono-dependencies-with-intel-mkl.sif cmake ../ -G "Ninja" \
        -DCMAKE_BUILD_TYPE=Release \
        -DENABLE_MODULE_POSTPROCESS=TRUE \
        -DENABLE_MODULE_PYTHON=TRUE \
        -DENABLE_MODULE_COSIMULATION=TRUE \
        -DENABLE_MODULE_IRRLICHT=TRUE \
        -DENABLE_MODULE_VEHICLE=FALSE \
        -DENABLE_MODULE_MULTICORE=TRUE \
        -DENABLE_MODULE_OPENGL=TRUE \
        -DENABLE_MODULE_FSI=TRUE \
        -DENABLE_MODULE_SYNCHRONO=TRUE \
        -DENABLE_MODULE_CSHARP=FALSE \
        -DENABLE_MODULE_GPU=TRUE \
        -DENABLE_MODULE_DISTRIBUTED=TRUE \
        -DENABLE_HDF5=TRUE \
        -DCMAKE_C_COMPILER=/usr/bin/gcc \
        -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
        -DCUDA_HOST_COMPILER=/usr/bin/gcc \
        -DPYTHON_EXECUTABLE=/usr/bin/python3 \
        -DEIGEN3_INCLUDE_DIR=/usr/include/eigen3 \
        -DCMAKE_VERBOSE_MAKEFILE=TRUE \
        -DChrono_DIR=/chrono-concrete/build/cmake \
        -DChrono_DATA_DIR=/chrono-concrete/build/data \
        -DCHRONO_CXX_FLAGS=TRUE\
        -DCHORONO_C_FLAGS=TRUE
pwd
# build project chrono-concrete
singularity exec --pwd /RobotGcode/build -B /projects/p31861/Users/LiuJiaqi/workdir/RobotGcode:/RobotGcode -B /projects/p31861/Users/AkpanWisdom/NEWCHRONO/chrono-concrete:/chrono-concrete /projects/p31861/Users/LaleErol/project-chrono-dependencies-with-intel-mkl.sif ninja -j 4

# Run job
cd /projects/p31861/Users/LiuJiaqi/workdir/RobotGcode

singularity exec -B /projects/p31861/Users/LiuJiaqi/workdir/RobotGcode:/RobotGcode -B /projects/p31861/Users/AkpanWisdom/NEWCHRONO/chrono-concrete:/chrono-concrete /projects/p31861/Users/LaleErol/project-chrono-dependencies-with-intel-mkl.sif /RobotGcode/build/RobotGCode
