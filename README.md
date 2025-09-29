# Howitt-Clower Economic Simulation

An object-oriented C++17 implementation of the Howitt-Clower economic simulation model, designed to study the emergence of monetary exchange in decentralized market economies.

## Introduction

This simulation implements the economic model described in the seminal work by Peter Howitt and Robert Clower, which explores how money can emerge spontaneously in a decentralized trading economy without central coordination. The model demonstrates how certain goods can become widely accepted as media of exchange (money) through the self-organizing behavior of individual traders seeking to maximize their utility.

## The Experiment

The Howitt-Clower model simulates an economy with multiple goods where agents:

- **Produce one good** and **desire a different good** (specialization)
- **Search for trading opportunities** in a decentralized marketplace
- **Open shops** as intermediaries when profitable
- **Adapt prices** based on realized income versus targets
- **Exit unprofitable ventures** stochastically

### Key Economic Mechanisms

1. **Search and Matching**: Agents sample a limited number of shops and adopt the best trading relationships
2. **Entry and Exit**: Entrepreneurs open new shops when research suggests profitability; unprofitable shops exit
3. **Price Formation**: Shop owners set prices based on adaptive income targeting
4. **Network Effects**: Agents discover trading opportunities through social connections (comrades and soulmates)

### Simulation Parameters

| Parameter | Value | Description |
|-----------|--------|-------------|
| Goods (n) | 10 | Number of different goods in the economy |
| Traders per type (bsize) | 24 | Number of agents for each (producer, consumer) pair |
| Total traders (m) | 2,160 | Total number of trading agents |
| Shop locations (K) | 200 | Maximum number of shops that can operate |
| Simulation length (T) | 20,000 | Number of weeks simulated |
| Search probability (λ) | 0.05 | Probability of searching when having trading links |
| Price adjustment (α) | 0.25 | Speed of adaptive price adjustment |
| Exit probability (θ) | 0.01 | Probability of unprofitable shop exit |
| Transaction cost (C) | 5.0 | Fixed cost per transaction |

## Original Research

This implementation is based on the research paper:

**Howitt, Peter, and Robert Clower. "The Emergence of Economic Organization." *Journal of Economic Theory* 94, no. 1 (2000): 1-31.**

The original paper introduces a theoretical framework for understanding how monetary institutions can arise endogenously in decentralized economies through individual optimization behavior, without requiring central authority or pre-existing monetary systems.

## Project Structure

```
howitt-clower/
├── src/                    # Source files (.cpp)
│   ├── main.cpp           # Program entry point
│   ├── PCG32.cpp          # Random number generator implementation
│   ├── Shop.cpp           # Shop class implementation
│   ├── Trader.cpp         # Trader class implementation
│   └── Simulation.cpp     # Main simulation logic
├── include/               # Header files (.h)
│   ├── PCG32.h            # RNG class definition
│   ├── ResearchResults.h  # Research outcome structure
│   ├── Shop.h             # Shop class definition
│   ├── Trader.h           # Trader class definition
│   └── Simulation.h       # Simulation class definition
├── dist/                  # Build output directory
│   ├── howclow_oop       # Executable
│   └── obj/              # Object files
├── original/              # Original reference implementations
├── Makefile              # Build configuration
└── README.md             # This file
```

## Dependencies

### System Requirements
- **C++ Compiler**: GCC 7.0+ or Clang 5.0+ with C++17 support
- **Make**: GNU Make or compatible build system
- **Operating System**: Linux, macOS, or Windows (with appropriate toolchain)

### Standard Library Dependencies
The project uses only standard C++ libraries:
- `<iostream>`, `<fstream>` - Input/output operations
- `<vector>`, `<array>` - Container classes
- `<algorithm>`, `<numeric>` - Standard algorithms
- `<functional>` - Function objects and utilities
- `<cstdint>`, `<cstdio>`, `<cstdlib>` - C standard library compatibility
- `<cmath>`, `<ctime>` - Mathematical and time functions
- `<string>`, `<cstring>` - String handling

**No external dependencies required** - the project is completely self-contained.

## Compilation and Usage

### Quick Start

1. **Clone the repository**:
   ```bash
   git clone https://github.com/deodoro/howitt-clower.git
   cd howitt-clower
   ```

2. **Compile the project**:
   ```bash
   make
   ```

3. **Run the simulation**:
   ```bash
   ./dist/howclow_oop
   ```

### Build Options

- **Clean build**: `make clean && make`
- **Clean only**: `make clean`
- **Debug build**: Modify `CXXFLAGS` in Makefile to include `-g -DDEBUG=1`

### Output Files

The simulation generates:
- **evol.fil**: Detailed simulation results and statistics
- **Console output**: Real-time progress indicators and summary statistics

## Main Features

### Core Components

1. **PCG32 Random Number Generator**
   - High-quality pseudorandom number generation
   - Reproducible results with configurable seeds
   - Uniform integer and floating-point distributions

2. **Shop Class**
   - Manages trading posts for goods exchange
   - Adaptive price setting based on income targets
   - Profit calculation and viability assessment
   - Support for both direct barter and indirect (monetary) exchange

3. **Trader Class**
   - Individual economic agents with production/consumption specialization
   - Social network connections (comrades producing same good, soulmates consuming same good)
   - Utility maximization through shop selection
   - Link management for trading relationships

4. **Simulation Engine**
   - Complete economic system orchestration
   - Weekly cycles of entry, matching, trading, and exit
   - Statistical analysis and monetary emergence detection
   - Configurable parameters for experimental variation

### Economic Analysis Features

- **Participation Rate Tracking**: Monitors fraction of agents successfully engaged in trade
- **Monetary Emergence Detection**: Identifies when specific goods become widely used as media of exchange
- **Price Convergence Analysis**: Measures deviation from theoretical equilibrium prices
- **Network Structure Evolution**: Tracks the development of trading relationships over time
- **Surplus Calculations**: Computes consumer and producer surplus measures

### Technical Features

- **Object-Oriented Design**: Clean separation of economic entities and behaviors
- **Memory Efficient**: Optimized data structures for large-scale simulation
- **Configurable Parameters**: Easy adjustment of economic and computational settings  
- **Comprehensive Logging**: Detailed output for research analysis
- **Cross-Platform**: Compatible with major operating systems

## Experimental Usage

### Parameter Modification

Key simulation parameters can be adjusted in `include/Simulation.h`:

```cpp
static constexpr int n = 10;           // Number of goods
static constexpr int bsize = 24;       // Traders per type
static constexpr int K = 200;          // Maximum shops
static constexpr int T = 20000;        // Simulation length
static constexpr double lambda = 0.05; // Search probability
static constexpr double alpha = 0.25;  // Price adjustment speed
static constexpr double theta = 0.01;  // Exit probability
```

### Research Applications

This simulation can be used to study:
- **Monetary emergence conditions**: Under what parameters do specific goods become money?
- **Market efficiency**: How does decentralized search affect economic welfare?
- **Network effects**: How do social connections influence market structure?
- **Policy interventions**: Effects of transaction costs, regulations, or subsidies

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

See the [GNU General Public License](https://www.gnu.org/licenses/) for more details.

## Author

**Jose Deodoro** <deodoro.filho@gmail.com> or <jdeoliv@gmu.edu>

## Acknowledgments

- Original theoretical framework by Peter Howitt and Robert Clower
- PCG random number generator algorithm by Melissa E. O'Neill
- Economic modeling insights from the agent-based computational economics literature
