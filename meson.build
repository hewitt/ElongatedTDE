project('ElongatedTDE.cpp', 'cpp')
# this is a Cpp project
compiler = meson.get_compiler('cpp')

add_global_arguments('-O2', language : 'cpp')
add_global_arguments('-DPETSC_Z', language : 'cpp')
#add_global_arguments('-DPARANOID', language : 'cpp')
#add_global_arguments('-DSLEPC', language : 'cpp')

lapack_dep = dependency('lapack', required:false)
blas_dep = dependency('blas', required:false)
petsc_dep = dependency('PETSc')
#slepc_dep = dependency('SLEPc')
#mpich_dep = dependency('mpich')
cppnoddy_dep = dependency('CppNoddy')

exe = executable('ElongatedTDE', 'ElongatedTDE.cpp', 
    dependencies : [ cppnoddy_dep, lapack_dep, blas_dep, petsc_dep ])

