import glob
import os

mode = ARGUMENTS.get('mode','gcc_release')

if mode=='gcc_release':
    compiler = 'g++'
    cflags = ' -Wfatal-errors -O3 '
    cflags += '-Werror -Wall -Wextra -pedantic -Wfatal-errors -Weffc++ -Wshadow -Wmissing-declarations -Wno-long-long -Wno-effc++ -Wno-parentheses -Wno-reorder -Wno-shadow -Wconversion'
elif mode=='gcc_debug':
    compiler = 'g++'
    cflags = ' -Wfatal-errors -O0 -g -pg '
elif mode=='clang_release':
    compiler= 'clang++'
    cflags = '-Weverything -Werror -ferror-limit=1 -Wno-error=padded -O3 '
else:
    raise NameError('unsupported mode')

source_dir = '.'
build_dir = 'build/'+mode
env = Environment(ENV = os.environ,
                  CXX=compiler,
                  CPPPATH=source_dir,
                  LIBPATH=['.',os.environ['HDF5_LIB_PATH']],
                  LIBS=['hdf5','hdf5_cpp'],
                  CXXFLAGS=cflags)
env.VariantDir(build_dir,source_dir)
env.Program(build_dir+'/tde',
            Glob(build_dir+'/*.cpp'))
