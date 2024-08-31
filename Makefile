.ONESHELL:

raytracer.o:
	cc -c -o raytracer.o -DRAVE_CUSTOM_VOXEL_TYPE=int src/raytracer.c -Ofast -msse4.2
test: tests/test.c dependencies/tgafunc.c raytracer.o
	cc -o test tests/test.c raytracer.o dependencies/tgafunc.c -Ofast -msse4.2 -lm -pthread
ifeq ($(OS),Windows_NT)
	.\test
else
	./test
endif

fractal: tests/fractal.c dependencies/tgafunc.c raytracer.o
	cc -o fractal tests/fractal.c dependencies/tgafunc.c -Ofast -msse4.2 -lm -pthread
ifeq ($(OS),Windows_NT)
	.\fractal
else
	./fractal
endif

clean: 
ifeq ($(OS),Windows_NT)
	del "*.o"  
	del "test"  
	del "fractal"  
else
	rm -R *.o
	rm -R test
	rm -R fractal
endif