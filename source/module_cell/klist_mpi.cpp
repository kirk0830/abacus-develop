/*
Note: on the plan of removal of GlobalC::parallel_kpoints
it is totally unreasonable to have such a class, it is very ridiculous that kpoints itself
cannot divide the MPI comm world.
Path: source/module_cell/parallel_kpoints.h, after refactor this class will be removed and
the functions will be moved to module_cell/klist_mpi.cpp

parallel_kpoints is instantiated in global.cpp, it's life span is presently unknown.
*/
