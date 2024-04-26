project('caltop_gm', ['c', 'fortran'], default_options : ['c_std=c11', 'cpp_std=c++14', 'fortran_std=f95'])

# Compiler and linker options
c_flags = ['-Wall', '-g', '-O2', '-fopenmp', '-DARCH="Linux"', '-DARPACK', '-DPARDISO', '-DMATRIXSTORAGE', '-DUSE_MT=1', '-DMKL_LP64']
f_flags = ['-Wall', '-g', '-O2', '-w', '-fno-second-underscore', '-fcray-pointer', '-x', 'f77-cpp-input', '-fopenmp']
link_args = [
    '-Wl,--start-group',
    '/opt/intel/oneapi/mkl/2022.2.0/lib/intel64/libmkl_gf_lp64.a',
    '/opt/intel/oneapi/mkl/2022.2.0/lib/intel64/libmkl_gnu_thread.a',
    '/opt/intel/oneapi/mkl/2022.2.0/lib/intel64/libmkl_core.a',
    '/opt/intel/oneapi/mkl/2022.2.0/lib/intel64/libmkl_blacs_openmpi_lp64.a',
    '-Wl,--end-group',
    '-lgomp', '-lpthread', '-lm', '-ldl'
]

# Includes and libraries
inc_dirs = include_directories('include')
mkl_include = include_directories('/opt/intel/oneapi/mkl/2022.2.0/include')
dependencies = [
    dependency('gfortran'),
    dependency('openmp')
]

# Source files
sources = files('src/main.c', 'src/other_c_files.c', 'src/other_fortran_files.f')

# Executable
executable('caltop_gm', sources, include_directories : [inc_dirs, mkl_include], c_args : c_flags, f_args : f_flags, link_args : link_args, dependencies : dependencies)

# Custom script to update date in source files (replacing the functionality of date.pl)
run_target('update_date',
    command : [find_program('date.pl')],
    depends : executable('caltop_gm')
)

# Cleaning up .old files if any
clean_files = find_program('clean_script.sh')
custom_target('clean_old_files',
    command : [clean_files],
    depend_files : [executable('caltop_gm')]
)
