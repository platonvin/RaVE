.ONESHELL:

test: tests/test.c
	cc -Ofast -msse4.2 tests/test.c dependencies/tgafunc.c src/raytracer_simd.c -pipe -o test
	timecmd test
test1: tests/test1.c
	cc -Ofast -march=native tests/test1.c -pipe -o test1
	test1
fractal: tests/fractal.c
	cc -Ofast -msse4.2 tests/fractal.c dependencies/tgafunc.c src/raytracer_simd.c -pipe -o fractal
	timecmd fractal