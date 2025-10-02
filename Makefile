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
TEST_RES_TARGET = $(DISTDIR)/test_research
TEST_MATCH_TARGET = $(DISTDIR)/test_match

# Files
TARGET = $(DISTDIR)/howclow_oop
SOURCES = $(filter-out $(SRCDIR)/TestResearch.cpp $(SRCDIR)/TestMatch.cpp, $(wildcard $(SRCDIR)/*.cpp))
TEST_RES_SOURCES = $(SRCDIR)/TestResearch.cpp
TEST_MATCH_SOURCES = $(SRCDIR)/TestMatch.cpp
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
TEST_RES_OBJECTS = $(TEST_RES_SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o) $(OBJDIR)/PCG32.o $(OBJDIR)/Shop.o $(OBJDIR)/Trader.o
TEST_MATCH_OBJECTS = $(TEST_MATCH_SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o) $(OBJDIR)/PCG32.o $(OBJDIR)/Shop.o $(OBJDIR)/Trader.o
HEADERS = $(wildcard $(INCDIR)/*.h)

.PHONY: all clean dirs test both

all: dirs $(TARGET)

both: dirs $(TARGET) $(TEST_RES_TARGET) $(TEST_MATCH_TARGET)

test: dirs $(TEST_RES_TARGET)

test-match: dirs $(TEST_MATCH_TARGET)

dirs:
	@mkdir -p $(DISTDIR) $(OBJDIR)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(TEST_RES_TARGET): $(TEST_RES_OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(TEST_MATCH_TARGET): $(TEST_MATCH_OBJECTS)
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
