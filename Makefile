.ONESHELL:

#example with object file
test: tests/test.c
	cc -Ofast -msse4.2 -DRAVE_CUSTOM_VOXEL_TYPE=int src/raytracer.c -c -o raytracer.o
	cc -Ofast -msse4.2 tests/test.c raytracer.o dependencies/tgafunc.c -o test
	test
#example with just .c inclusion
fractal: tests/fractal.c
	cc -Ofast -msse4.2 tests/fractal.c dependencies/tgafunc.c -o fractal
	fractal