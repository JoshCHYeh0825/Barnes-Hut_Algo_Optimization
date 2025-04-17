CC = gcc
CFLAGS = -Wall -Wextra -O2
LDFLAGS = -lgdi32 -mwindows

SRCS = vec2.c body.c quadtree.c
OBJS = $(SRCS:.c=.o)
MAIN_OBJ = main.o
TARGET = barnes_hut.exe

all: $(TARGET)

$(TARGET): $(OBJS) $(MAIN_OBJ)
	$(CC) $(MAIN_OBJ) $(OBJS) -o $(TARGET) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	del /F /Q $(OBJS) $(MAIN_OBJ) $(TARGET)

run: $(TARGET)
	.\$(TARGET)

.PHONY: all clean run 