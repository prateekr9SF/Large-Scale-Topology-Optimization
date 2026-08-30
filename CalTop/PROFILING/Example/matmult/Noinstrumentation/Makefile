# ==============================================================
# Matrix Multiplication Example
# ==============================================================

# Compiler
CC = g++-12

# Compiler flags
CFLAGS = -O2 -g -Wall -Wextra

# Directories
SRC_DIR = src
INC_DIR = include
OBJ_DIR = obj
BIN_DIR = bin

# Include path
INCLUDES = -I$(INC_DIR)

# Source files
SOURCES = $(SRC_DIR)/main.c \
          $(SRC_DIR)/matrix.c

# Object files
OBJECTS = $(OBJ_DIR)/main.o \
          $(OBJ_DIR)/matrix.o

# Executable
TARGET = $(BIN_DIR)/matrix_multiply.exe


# ==============================================================
# Default target
# ==============================================================

.PHONY: all
all: $(TARGET)


# ==============================================================
# Link executable
# ==============================================================

$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(CC) $(OBJECTS) -o $(TARGET) -lm


# ==============================================================
# Compile source files
# ==============================================================

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@


# ==============================================================
# Create build directories
# ==============================================================

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)


# ==============================================================
# Run executable
# ==============================================================

.PHONY: run
run: $(TARGET)
	./$(TARGET)


# ==============================================================
# Clean build files
# ==============================================================

.PHONY: clean
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)


# ==============================================================
# Display compiler information
# ==============================================================

.PHONY: compiler-info
compiler-info:
	@echo "========================================"
	@echo "Compiler configuration"
	@echo "========================================"
	@echo "CC = $(CC)"
	@$(CC) --version | head -n 1
	@echo "========================================"