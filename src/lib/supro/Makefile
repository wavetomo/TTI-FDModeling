BUILD_DIR = ./build
TARGET = $(INSTALL_LIB_DIR)/libsupro.a
SRC = $(wildcard *.c)
OBJECTS = $(patsubst %.c, $(BUILD_DIR)/%.o, $(SRC))
CFLAGS = -Wall -std=c11 -qopenmp -fPIC -O3

$(TARGET): $(OBJECTS)
	ar -rv $@ $^
$(BUILD_DIR)/%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(TARGET)