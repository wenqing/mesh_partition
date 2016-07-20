# partmesh

==============
A tool for preparing the partitioned mesh for OpenGeoSys (OGS) 5 (https://github.com/ufz/ogs5) for parallel computing. METIS is used as the partitioning application. Once the source code is compiled, one gets two executable files, partmesh and mpmetis.

## Usage:
Assuming that a mesh save in file foo.msh is going to be partitioned, a partitioning can be conducted as follows for two different OGS 5 compilations for parallel FEM computing.

###### Compiled with Cmake option ``-DOGS_CONFIG=MPI`` or ``-DOGS_FEM_MPI=ON`` for old version
1. Convert the OGS mesh data to METIS input data
  
      ``partmesh --ogs2metis foo``
2. Partition mesh data

      ``[full path]/partmesh --metis2ogs -np 5 -s -e foo``

    where option ``-s`` means to call mpmetis via sysem call to partition the mesh, and option ``-e`` means element wise partitioning. Alternatively, mpmetis can be run beforehand to save memory usage for large mesh, and then run parmesh without -s option. For example,

      ``[full path]/mpmetis foo.mesh 5``

      ``[full path]/partmesh --metis2ogs -np 5 -e foo``
      
      Option ``-ogsmsh`` can be used here to output the partitioned mesh data into files for Gmsh for visualization.
    
    Note, in this step the command must be run with its full path.

###### Compiled with Cmake option ``-DOGS_CONFIG=PETSC`` or ``-DOGS_FEM_PETSC=ON`` for old version

1. Convert the OGS mesh data to METIS input data
  
    For modelling with linear elements (processes excluding deformation), run command

      ``partmesh --ogs2metis foo`` 

    For modelling with quadratic elements (processes including deformation, e.g. M, TM, HM, H2M), run command
      ``partmesh --ogs2metis -q foo`` 
   
    Option ``-q`` indicates quadratic element.
2. Partition mesh data

      ``[full path]/partmesh --metis2ogs -np 5 -s -n foo``

    where option ``-s`` means to call mpmetis via sysem call to partition the mesh, , and option ``-n`` means node wise partitioning. Option '-q' has to be selected if it is already used in the first step in converting OGS mesh data to METIS input data.
    
    Alternatively, mpmetis can be run beforehand to save memory usage for large mesh, and then run parmesh without -s option. For example,

      ``[full path]/mpmetis foo.mesh 5``

      ``[full path]/partmesh --metis2ogs -np 5 -n foo``
      
      Option ``-odom`` can be used here to output the partitioned mesh data into VTK files for visualization.
    
    Note, in this step the command must be run with its full path.

Always, the full options of the tool can be displayed by run command

      ``partmesh --help`` 
