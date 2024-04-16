project('mixed_language_project', ['c', 'fortran'], default_options : ['c_std=c11', 'fortran_std=f95'])

# Define executables
main_sources = ['main.c']  # Add C sources
fortran_sources = ['stop.f', 'closefile.f']  # Add Fortran sources

# Include directories
inc_dirs = include_directories('.')

# Dependencies
# If there are any specific libraries or dependencies, list them here
# e.g., deps = dependency('some_library')

# Create an executable
executable('my_program', sources : [main_sources, fortran_sources], include_directories : inc_dirs, link_with : [], dependencies : [])

