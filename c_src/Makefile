# Name of final executable
TARGET = barnes_hut

# Compiler settings
CC = gcc
NVCC = nvcc
CFLAGS = -O3 -Wall -I. -fopenmp
NVCCFLAGS = -O3 -I. -Xcompiler -fopenmp

# List of source files
CPU_SRCS = main.c body.c vec2.c quadtree.c
CUDA_SRCS = barnes_hut_CUDA.cu

# Object files
OBJS = $(CPU_SRCS:.c=.o) $(CUDA_SRCS:.cu=.o)

# Libraries
LIBS = -lm -lSDL

# Default rule
all: $(TARGET)

# Rule to build final executable
$(TARGET): $(OBJS)
	$(NVCC) $(NVCCFLAGS) -o $@ $(OBJS) $(LIBS)

# Rule for C files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Rule for CUDA files
%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(TARGET) $(OBJS)
