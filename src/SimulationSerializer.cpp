#include "SimulationSerializer.h"
#include <iomanip>

// Helper to serialize SimulationInfo to JSON
static void serializeSimulationInfo(const SimulationInfo& info, std::ostream& out) {
    out << info.to_json();
}

// Helper to serialize RunInfo to JSON
static void serializeRunInfo(const RunInfo* run, std::ostream& out) {
    out << run->to_json();
}

void SimulationSerializer::serialize(const SimulationInfo& info, const std::vector<std::vector<RunInfo*>>& runs, std::ostream& out) {
    out << "{\n";
    out << "  \"SimulationInfo\": ";
    serializeSimulationInfo(info, out);
    out << ",\n  \"runs_per_slope\": [\n";
    for (size_t i = 0; i < runs.size(); ++i) {
        out << "    [\n";
        for (size_t j = 0; j < runs[i].size(); ++j) {
            serializeRunInfo(runs[i][j], out);
            if (j + 1 < runs[i].size()) out << ",\n";
        }
        out << "    ]";
        if (i + 1 < runs.size()) out << ",\n";
    }
    out << "\n  ]\n";
    out << "}";
}
