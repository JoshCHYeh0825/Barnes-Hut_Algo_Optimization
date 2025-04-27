CC = gcc
CFLAGS = -Wall -Wextra -O2 -fopenmp `sdl-config --cflags`
LDFLAGS = -fopenmp `sdl-config --libs` -lm

SRCS = vec2.c body.c quadtree.c
OBJS = $(SRCS:.c=.o)
MAIN_OBJ = main.o
TARGET = barnes_hut

all: $(TARGET)

$(TARGET): $(OBJS) $(MAIN_OBJ)
	$(CC) $(MAIN_OBJ) $(OBJS) -o $(TARGET) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(MAIN_OBJ) $(TARGET)

run: $(TARGET)
	./$(TARGET)

.PHONY: all clean run
