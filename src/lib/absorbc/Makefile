BUILD_DIR = ./build

TARGET = $(INSTALL_LIB_DIR)/libabsorbc.a
SRC = $(wildcard *.c)
OBJECTS = $(patsubst %.c, $(BUILD_DIR)/%.o, $(SRC))

$(TARGET): $(OBJECTS)
	ar -rv $@ $^
$(BUILD_DIR)/%.o: %.c
	$(CC) $(LIBFLAG) $(INC) -c $< -o $@

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(TARGET)