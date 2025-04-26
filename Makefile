CC = gcc
CFLAGS = -Wall -Wextra -O2 -I/opt/homebrew/include
LDFLAGS = -L/opt/homebrew/lib -lSDL2 -lm
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