.ONESHELL:

test: tests/test.c
	cc -Ofast -msse4.2 tests/test.c dependencies/tgafunc.c src/raytracer.c -pipe -o test
	timecmd test
test1: tests/test1.c
	cc -Ofast -march=native -lpthread tests/test1.c -pipe -o test1
	test1
fractal: tests/fractal.c
	cc -Ofast -msse4.2 -lpthread tests/fractal.c dependencies/tgafunc.c src/raytracer.c -pipe -o fractal
	timecmd fractal