CC = gcc

CFLAGS = -Wall -Wextra -O3 -ffast-math `sdl-config --cflags` -fopenmp
LDFLAGS = `sdl-config --libs` -lm -fopenmp

SRCS = vec2.c body.c quadtree.c
OBJS = $(SRCS:.c=.o)
MAIN_OBJ = main.o
TARGET = barnes_hut

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(OBJS) $(MAIN_OBJ)
	$(CC) $(MAIN_OBJ) $(OBJS) -o $(TARGET) $(LDFLAGS)
	rm -f $(OBJS) $(MAIN_OBJ)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(MAIN_OBJ) $(TARGET)

run: $(TARGET)
	./$(TARGET)