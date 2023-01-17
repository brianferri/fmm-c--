CXX := g++

ifeq ($(OS),Windows_NT)
	detected_OS := Windows
	RM := rmdir /Q /S
	UN := del
	FixPath = $(subst /,\,$1)
	ext := .exe
else
	detected_OS := $(shell uname)
	RM := rm -f -r
	UN := $(RM)
	FLAGS := -p
	FixPath = $1
endif

MD := mkdir $(FLAGS)

# Dependency flags
CXXFLAGS += -MMD -MP

# Compiler flags
CXXFLAGS += -std=c++17 -pedantic-errors
CXXFLAGS += -Wall -Wextra
CXXFLAGS += -g3 -O3

# Linker flags
CPPFLAGS += -Iinclude

DIRS := $(patsubst src/%, %, $(wildcard src/*))
PROG_SOURCES := $(wildcard src/*/*.cpp)
OBJECTS := $(patsubst src/%.cpp, build/%.o, $(PROG_SOURCES))
DEPENDENCIES := $(patsubst src/%.cpp, build/%.d, $(PROG_SOURCES))

all: prepare dirs main run clean
	
-include $(call FixPath,$(DEPENDENCIES))

prepare:
	@echo Preparing build directories
	-@$(RM) build

dirs:
	@$(MD) $(call FixPath,$(patsubst %, build/%, $(DIRS)))
	@echo Created directories: $(DIRS)

build/%.o: src/%.cpp
	@echo Compiling: $< : $@
	@$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(call FixPath,$<) -o $(call FixPath,$@)

main: $(OBJECTS)
	@echo Linking $@
	@$(CXX) $(CXXFLAGS) -o $@ $(call FixPath,$^) $(LDFLAGS)

run: main$(ext)
	@echo Running main
	@./main

clean:
	@echo Cleaning up executable
	-@$(UN) $(call FixPath,./main$(ext))

.PHONY: show dirs