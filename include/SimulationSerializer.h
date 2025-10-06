#ifndef SIMULATION_SERIALIZER_H
#define SIMULATION_SERIALIZER_H

#include "Simulation.h"
#include "SimulationInfo.h"
#include "RunInfo.h"
#include <ostream>
#include <vector>

class SimulationSerializer {
public:
    /**
     * Serializes the Simulation's SimulationInfo and runInfo vectors to a JSON stream.
     * @param sim The Simulation to serialize.
     * @param out The output stream to write JSON to.
     */
    static void serialize(const SimulationInfo& info, const std::vector<std::vector<RunInfo*>>& runs, std::ostream& out);
};

#endif // SIMULATION_SERIALIZER_H
