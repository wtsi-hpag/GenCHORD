## This is a custom makefile built with the help of an LLM because I don't know how makefiles work and have been using the same one for 10 years. 
## I am told that it robustly compiles two distinct executables (MAIN and PREPROCESS) at the same time. 
## Fingers crossed.

# Project Names
MAIN_PROJECT = genchord
PREPROCESS_PROJECT = depthsplitter

# Compiler
CC = g++

# Compiler and Linker Options
CXXFLAGS = -std=c++20 -pthread -O3 -w -Ilibs/JSL
LDFLAGS = -lpthread

# Dependency Flags
DEPFLAGS = -MMD -MP

# Source and Build Directories
MAIN_SRC_DIR = src
PREPROCESS_SRC_DIR = samdepth
BUILD_DIR = build

### DON'T EDIT BELOW HERE


# Find all source files
MAIN_SOURCE_FILES := $(shell find $(MAIN_SRC_DIR) -name "*.cpp")
PREPROCESS_SOURCE_FILES := $(shell find $(PREPROCESS_SRC_DIR) -name "*.cpp")

# Generate object files
MAIN_OBJECTS := $(patsubst $(MAIN_SRC_DIR)/%.cpp, $(BUILD_DIR)/main/%.o, $(MAIN_SOURCE_FILES))
PREPROCESS_OBJECTS := $(patsubst $(PREPROCESS_SRC_DIR)/%.cpp, $(BUILD_DIR)/preprocess/%.o, $(PREPROCESS_SOURCE_FILES))

# Dependency Files
MAIN_DEPS := $(MAIN_OBJECTS:.o=.d)
PREPROCESS_DEPS := $(PREPROCESS_OBJECTS:.o=.d)

# Default Target: Build both projects
.PHONY: all
all: $(MAIN_PROJECT) $(PREPROCESS_PROJECT)

# Build Main Project
$(MAIN_PROJECT): $(MAIN_OBJECTS)
	@echo "Linking $(MAIN_PROJECT)..."
	$(CC) -o $@ $^ $(LDFLAGS)

# Build Preprocess Project
$(PREPROCESS_PROJECT): $(PREPROCESS_OBJECTS)
	@echo "Linking $(PREPROCESS_PROJECT)..."
	$(CC) -o $@ $^ $(LDFLAGS)

# Compile Main Source Files
$(BUILD_DIR)/main/%.o: $(MAIN_SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	@echo "Compiling $< for $(MAIN_PROJECT)..."
	$(CC) $(CXXFLAGS) $(DEPFLAGS) -c -o $@ $<

# Compile Preprocess Source Files
$(BUILD_DIR)/preprocess/%.o: $(PREPROCESS_SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	@echo "Compiling $< for $(PREPROCESS_PROJECT)..."
	$(CC) $(CXXFLAGS) $(DEPFLAGS) -c -o $@ $<

# Include Dependencies
-include $(MAIN_DEPS)
-include $(PREPROCESS_DEPS)

# Run Main Project
.PHONY: run-main
run-main: $(MAIN_PROJECT)
	./$(MAIN_PROJECT)

# Run Preprocess Project
.PHONY: run-preprocess
run-preprocess: $(PREPROCESS_PROJECT)
	./$(PREPROCESS_PROJECT)

# Clean Targets
.PHONY: clean
clean:
	@echo "Cleaning up..."
	rm -rf $(MAIN_PROJECT) $(PREPROCESS_PROJECT) $(BUILD_DIR)

.PHONY: depclean
depclean:
	@echo "Removing dependency files..."
	rm -f $(MAIN_DEPS) $(PREPROCESS_DEPS)

.PHONY: clean-all
clean-all: clean depclean
