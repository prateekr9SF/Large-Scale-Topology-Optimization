# ==============================================================
# Matrix Multiplication Example
# ==============================================================
TAU_ROOT = /home/prateek/Documents/Software/TAU/TAU_install
# Compiler
CC = gcc-12
TAU_CC = tau_cc.sh

# Compiler flags
CFLAGS = -O2 -g -Wall -Wextra -fopenmp -DPROFILING_ON
LDFLAGS = -fopenmp -lm

# Directories
SRC_DIR = src
INC_DIR = include
OBJ_DIR = obj
BIN_DIR = bin

# Include paths
INCLUDES = -I$(INC_DIR) \
           -I$(TAU_ROOT)/include

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
# Link executable with TAU
# ==============================================================

$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(TAU_CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS)


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
	rm -f profile.*


# ==============================================================
# Display compiler information
# ==============================================================

.PHONY: compiler-info
compiler-info:
	@echo "========================================"
	@echo "Compiler configuration"
	@echo "========================================"
	@echo "CC     = $(CC)"
	@$(CC) --version | head -n 1
	@echo
	@echo "TAU_CC = $(TAU_CC)"
	@$(TAU_CC) -show
	@echo "========================================"