BUILD_DIR = ./build
SRC = $(wildcard *.c)
TARGET = $(patsubst %.c, $(INSTALL_BIN_DIR)/%, $(SRC))
OBJECTS = $(patsubst %.c, $(BUILD_DIR)/%.o, $(SRC))

LFLAGS = -lfdmod -labsorbc -lrockmod -lsupro \
-lsu -lpar -lcwp -lm

all:$(TARGET)
$(TARGET): $(OBJECTS)
	$(CC) $(EXEFLAG) \
	$(patsubst %, $(BUILD_DIR)/%.o, $(notdir $@)) \
	$(LIB) $(LFLAGS) -o $@
$(BUILD_DIR)/%.o: %.c
	$(CC) $(EXEFLAG) $(INC) -c $< -o $@

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(TARGET)