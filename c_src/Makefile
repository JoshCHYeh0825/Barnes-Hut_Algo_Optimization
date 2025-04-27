# Compiler
NVCC = nvcc
CFLAGS = -O3 -Xcompiler -fopenmp

# Libraries
LIBS = -lSDL

# Executable
TARGET = barnes_hut_sim

# Objects to compile
OBJECTS = main.o body.o quadtree.o vec2_host.o barnes_hut_CUDA.o

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(NVCC) $(CFLAGS) -o $@ $^ $(LIBS)

%.o: %.c
	$(NVCC) $(CFLAGS) -c $< -o $@

%.o: %.cu
	$(NVCC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGET)
