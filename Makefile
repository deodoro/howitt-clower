CXX = g++
CXXFLAGS = -std=gnu++17 -O2 -pipe -Wall -Wextra
TARGET = howclow_oop
SOURCES = main.cpp PCG32.cpp Shop.cpp Trader.cpp Simulation.cpp
OBJECTS = $(SOURCES:.cpp=.o)
HEADERS = PCG32.h Shop.h Trader.h Simulation.h ResearchResults.h

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

# Dependencies
main.o: main.cpp Simulation.h
PCG32.o: PCG32.cpp PCG32.h
Shop.o: Shop.cpp Shop.h
Trader.o: Trader.cpp Trader.h Shop.h PCG32.h
Simulation.o: Simulation.cpp Simulation.h PCG32.h Shop.h Trader.h ResearchResults.h
