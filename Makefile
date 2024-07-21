.ONESHELL:

debug_Flags = -pipe

# TODO -findirect-inlining
test: tests/test.c
	cc -Ofast -fassociative-math -fopt-info-vec-missed -msse4 -findirect-inlining tests/test.c dependencies/tgafunc.c src/raytracer.c -o test
	timecmd test
