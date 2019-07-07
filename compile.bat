gfortran -c m_general.f90 -g
gfortran -c m_shell.f90 -g
gfortran -c m_link.f90 -g
gfortran -c main.f90 -g
gfortran -o TermSymbolProgram *.o -g