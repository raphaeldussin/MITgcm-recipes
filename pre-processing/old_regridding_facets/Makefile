
all: mod_akima_1d.so mod_drown_sosie.so

mod_akima_1d.so: mod_akima_1d.f90
	f2py -c -m mod_akima_1d mod_akima_1d.f90

mod_drown_sosie.so: mod_drown_sosie.f90
	f2py -c -m mod_drown_sosie mod_drown_sosie.f90

clean:
	\rm -Rf *.so *.dSYM
