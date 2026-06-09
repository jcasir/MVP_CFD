# Makefile for the MVP_CFD solver 

CXX      = g++
STANDARD = -std=c++20
WARN     = -Wall -Wextra -Wpedantic

# ── Directory ──────────────────────────────────────────────────────────────── 
SOLDIR   = solver
SRCDIR   = src

# ── Sources ──────────────────────────────────────────────────────────────────
SOURCES  = $(wildcard $(SRCDIR)/$(SOLDIR)/*.cpp) $(wildcard $(SRCDIR)/*.cpp)

# ── Default configuration (default: release) ───────────────────────────────── 
BUILD    ?= release

ifeq ($(BUILD),debug)
  OBJDIR    = obj/debug
  CXXFLAGS  = $(STANDARD) $(WARN) -g -O0 -Iinclude
  EXECSUFFIX = _debug
else ifeq ($(BUILD),asan)
  OBJDIR    = obj/asan
  CXXFLAGS  = $(STANDARD) $(WARN) -g -O1 -fsanitize=address,undefined -Iinclude
  LDFLAGS   = -fsanitize=address,undefined
  EXECSUFFIX = _asan
else
  OBJDIR    = obj/release
  CXXFLAGS  = $(STANDARD) $(WARN) -O3 -march=native -Iinclude
  EXECSUFFIX =
endif

EXECUTABLE = solver$(EXECSUFFIX)
CONFIG     = config.cfg

# ── Objects and automatic dependencies ─────────────────────────────────────── 
OBJECTS = $(patsubst %.cpp, $(OBJDIR)/%.o, $(SOURCES))
DEPFILES = $(OBJECTS:.o=.d)

# ── Main Target ──────────────────────────────────────────────────────────────
.PHONY: all run clean cleanall debug asan help

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $^ -o $@ $(LDFLAGS)
	@echo "✓ Build [$(BUILD)] completed → $(EXECUTABLE)"

# ── Compilation with automatic dependencies (-MMD -MP) ─────────────────────── 
$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(DEPFILES)

# ── Shortcuts for alternative builds ─────────────────────────────────────────
debug:
	$(MAKE) BUILD=debug

asan:
	$(MAKE) BUILD=asan

# ── Execution ────────────────────────────────────────────────────────────────
run: all
	./$(EXECUTABLE) $(CONFIG)

# ── Cleaning ───────────────────────────────────────────────────────────────── 
clean:
	rm -rf obj/
	rm -f solver solver_debug solver_asan
	rm -f results/*.csv
	rm -f results/*.vtu
	rm -f results/*.pvd

cleanall: clean
	rm -f *.dat *.txt

# ── Help ────────────────────────────────────────────────────────────────────── 
help:
	@echo ""
	@echo "  Makefile — MVP_CFD Solver"
	@echo ""
	@echo "  make [BUILD=release]   Compile in optimized mode (default)"
	@echo "  make BUILD=debug       Compile with debug symbols"
	@echo "  make BUILD=asan        Compile with AddressSanitizer" 
	@echo "  make run               Compile and run" 
	@echo "  make clean             Removes all objects and executables"
	@echo "  make cleanall          Also removes .dat/.txt files" 
	@echo "  make help              Show this message" 
	@echo ""