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
  EXECSUFFIX = _asan [cite: 2]
else
  OBJDIR    = obj/release
  CXXFLAGS  = $(STANDARD) $(WARN) -O3 -march=native -Iinclude
  EXECSUFFIX =
endif

EXECUTABLE = solver1d$(EXECSUFFIX)
CONFIG     = config.cfg

# ── Objects and automatic dependencies ─────────────────────────────────────── 
OBJECTS  = $(addprefix $(OBJDIR)/, $(SOURCES:.cpp=.o))
DEPFILES = $(OBJECTS:.o=.d)

# ── Main Target ──────────────────────────────────────────────────────────────
.PHONY: all run clean cleanall plot debug asan help

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
	rm -f solver1d solver1d_debug solver1d_asan
	rm -f results/*.csv

cleanall: clean
	rm -f *.dat *.txt

# ── Plot with automatically detected files ─────────────────────────────────── 
DATFILES = $(wildcard *.dat)
plot:
ifeq ($(DATFILES),)
	@echo "No .dat files found. Run 'make run' first."
else
	@echo "$(foreach f,$(DATFILES),'$(f)' w l lw 2 title '$(f)', )" | \
	  gnuplot -e "set terminal qt; \
	              set grid; set xlabel 'x'; set ylabel 'u'; \
	              plot $(shell echo "$(foreach f,$(DATFILES),'$(f)' w l lw 2 title '$(f)',)" | sed 's/,$$//' ); \
	              pause -1"
endif

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
	@echo "  make plot              Visualize all .dat files with gnuplot" 
	@echo "  make help              Show this message" 
	@echo ""