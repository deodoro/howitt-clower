# Compiler and flags
CXX = g++
CXXFLAGS = -std=gnu++17 -O2 -pipe -Wall -Wextra -Iinclude

# Directories
SRCDIR = src
INCDIR = include
DISTDIR = dist
OBJDIR = $(DISTDIR)/obj
CXXFLAGS += -g -O0
CFLAGS   += -g -O0
LDFLAGS  += 

# Files
TARGET = $(DISTDIR)/howclow_oop
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
HEADERS = $(wildcard $(INCDIR)/*.h)

.PHONY: all clean dirs

all: dirs $(TARGET)

dirs:
	@mkdir -p $(DISTDIR) $(OBJDIR)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(DISTDIR)

# Specific dependencies
$(OBJDIR)/main.o: $(SRCDIR)/main.cpp $(INCDIR)/Simulation.h
$(OBJDIR)/PCG32.o: $(SRCDIR)/PCG32.cpp $(INCDIR)/PCG32.h
$(OBJDIR)/Shop.o: $(SRCDIR)/Shop.cpp $(INCDIR)/Shop.h
$(OBJDIR)/Trader.o: $(SRCDIR)/Trader.cpp $(INCDIR)/Trader.h $(INCDIR)/Shop.h $(INCDIR)/PCG32.h
$(OBJDIR)/Simulation.o: $(SRCDIR)/Simulation.cpp $(INCDIR)/Simulation.h $(INCDIR)/PCG32.h $(INCDIR)/Shop.h $(INCDIR)/Trader.h $(INCDIR)/ResearchResults.h
