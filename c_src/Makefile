# Compiler
NVCC = nvcc
CFLAGS = -O3 -Xcompiler -fopenmp

# Libraries
LIBS = -lSDL

# Executable
TARGET = barnes_hut_sim

# Sources
CPU_SRCS = main.c vec2.c body.c quadtree.c
GPU_SRCS = barnes_hut_CUDA.cu

# Object files
CPU_OBJS = $(CPU_SRCS:.c=.o)
GPU_OBJS = $(GPU_SRCS:.cu=.o)

# Default target
all: $(TARGET)

$(TARGET): $(CPU_OBJS) $(GPU_OBJS)
	$(NVCC) $(CFLAGS) -o $@ $^ $(LIBS)

%.o: %.c
	$(NVCC) $(CFLAGS) -c $< -o $@

%.o: %.cu
	$(NVCC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(CPU_OBJS) $(GPU_OBJS) $(TARGET)
