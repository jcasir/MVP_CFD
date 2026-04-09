# Makefile per il risolutore CFD 1D

CXX = g++
CXXFLAGS = -std=c++17 -Wall -O3 -march=native
LDFLAGS = 

# File sorgente
SOURCES = main.cpp AdvectionDiffusionSolver.cpp ConfigParser.cpp
HEADERS = AdvectionDiffusionSolver.h ConfigParser.h InitialConditions.h math_settings.h
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = solver1d

# Target principale
all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS)
	@echo "Compilazione completata! Eseguibile: $(EXECUTABLE)"

# Regola per compilare i file oggetto
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Esegui il programma
run: $(EXECUTABLE)
	./$(EXECUTABLE)

# Pulizia file compilati
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
	rm -f *.dat

# Pulizia completa (include anche i file di output)
cleanall: clean
	rm -f *.dat *.txt

# Visualizza i risultati con gnuplot (se disponibile)
plot:
	@echo "set terminal qt; \
	       set grid; \
	       set xlabel 'x'; \
	       set ylabel 'u'; \
	       plot 'avvezione_t0.dat' w l lw 2 title 'Avvezione t=0', \
	            'avvezione_t05.dat' w l lw 2 title 'Avvezione t=0.5'; \
	       pause -1" | gnuplot

# Debug build
debug: CXXFLAGS = -std=c++17 -Wall -g -O0
debug: clean $(EXECUTABLE)

# Help
help:
	@echo "Makefile per Solver CFD 1D"
	@echo ""
	@echo "Target disponibili:"
	@echo "  make          - Compila il programma"
	@echo "  make run      - Compila ed esegue"
	@echo "  make clean    - Rimuove file oggetto ed eseguibile"
	@echo "  make cleanall - Rimuove tutto (inclusi file .dat)"
	@echo "  make plot     - Visualizza risultati con gnuplot"
	@echo "  make debug    - Compila versione debug"
	@echo "  make help     - Mostra questo messaggio"

.PHONY: all run clean cleanall plot debug help
