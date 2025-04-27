# Compiler
NVCC = nvcc
CFLAGS = -O3 -Xcompiler -fopenmp

# Libraries
LIBS = -lSDL

# Executable
TARGET = barnes_hut_sim

CPU_OBJS = main.o body.o quadtree.o vec2_host.o
GPU_OBJS = barnes_hut_CUDA.o vec2_device.o

all: $(TARGET)

$(TARGET): $(CPU_OBJS) $(GPU_OBJS)
	$(NVCC) $(CFLAGS) -o $@ $^ $(LIBS)

%.o: %.c
	$(NVCC) $(CFLAGS) -c $< -o $@

%.o: %.cu
	$(NVCC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o $(TARGET)
