# Compiler
NVCC = nvcc
CFLAGS = -O3 -Xcompiler -fopenmp

# Libraries
LIBS = -lSDL

# Executable
TARGET = barnes_hut

# Object files
OBJECTS = main.o body.o quadtree.o vec2.o vec2_CUDA.o barnes_hut_CUDA.o

# Default target
all: $(TARGET)

# Link the final executable
$(TARGET): $(OBJECTS)
	$(NVCC) $(CFLAGS) -o $@ $^ $(LIBS)

# Compile .c files into .o files
%.o: %.c
	$(NVCC) $(CFLAGS) -c $< -o $@

# Compile .cu files into .o files
%.o: %.cu
	$(NVCC) $(CFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJECTS) $(TARGET)
