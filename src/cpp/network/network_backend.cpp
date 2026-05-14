#include "network_backend.h"
#include "haplotype.h"
#include "tcs_algorithm.h"
#include "msn_algorithm.h"
#include "mjn_algorithm.h"
#include "../third_party/json.hpp"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>

using json = nlohmann::json;

namespace haplokit {
namespace network {

std::string NetworkBackend::parse_algorithm(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <algorithm>\n";
        std::cerr << "Algorithms: tcs, msn, mjn\n";
        return "";
    }

    std::string algo = argv[1];
    std::transform(algo.begin(), algo.end(), algo.begin(), ::tolower);

    if (algo != "tcs" && algo != "msn" && algo != "mjn") {
        std::cerr << "Unknown algorithm: " << algo << "\n";
        std::cerr << "Supported: tcs, msn, mjn\n";
        return "";
    }

    return algo;
}

std::string NetworkBackend::read_stdin() {
    std::ostringstream buffer;
    buffer << std::cin.rdbuf();
    return buffer.str();
}

// Simple JSON parser for our specific format
NetworkBackend::InputData NetworkBackend::parse_input_json(const std::string& json_str) {
    InputData data;

    try {
        json j = json::parse(json_str);

        // Parse sequences array
        if (j.contains("sequences") && j["sequences"].is_array()) {
            for (const auto& seq : j["sequences"]) {
                if (seq.is_string()) {
                    data.sequences.push_back(seq.get<std::string>());
                }
            }
        }

        // Parse samples array (array of arrays)
        if (j.contains("samples") && j["samples"].is_array()) {
            for (const auto& sample_list : j["samples"]) {
                std::vector<std::string> samples;
                if (sample_list.is_array()) {
                    for (const auto& sample : sample_list) {
                        if (sample.is_string()) {
                            samples.push_back(sample.get<std::string>());
                        }
                    }
                }
                data.samples.push_back(samples);
            }
        }

        // Parse epsilon
        if (j.contains("epsilon")) {
            if (j["epsilon"].is_number_integer()) {
                data.epsilon = j["epsilon"].get<int>();
            }
        }

    } catch (const json::exception& e) {
        std::cerr << "JSON parse error: " << e.what() << "\n";
    }

    return data;
}

std::string NetworkBackend::generate_output_json(const OutputData& data) {
    json output = json::object();

    // Build nodes array
    json nodes = json::array();
    for (size_t i = 0; i < data.node_ids.size(); ++i) {
        json node = {
            {"id", data.node_ids[i]},
            {"sequence", data.node_sequences[i]},
            {"is_median", data.node_is_median[i]},
            {"samples", data.node_samples[i]}
        };
        nodes.push_back(node);
    }
    output["nodes"] = nodes;

    // Build edges array
    json edges = json::array();
    for (size_t i = 0; i < data.edge_sources.size(); ++i) {
        json edge = {
            {"source", data.edge_sources[i]},
            {"target", data.edge_targets[i]},
            {"weight", data.edge_weights[i]}
        };
        edges.push_back(edge);
    }
    output["edges"] = edges;

    return output.dump(2);  // Pretty print with 2-space indent
}

NetworkBackend::OutputData NetworkBackend::compute_network(
    const std::string& algorithm,
    const InputData& input) {

    OutputData output;

    // ES.20 + ES.23: Prepare haplotypes with uniform initialization
    std::vector<Haplotype> haplotypes;
    haplotypes.reserve(input.sequences.size());

    for (size_t i{0}; i < input.sequences.size(); ++i) {
        Haplotype hap{};
        hap.sequence = input.sequences[i];
        if (i < input.samples.size()) {
            hap.samples = input.samples[i];
        }
        haplotypes.push_back(std::move(hap));
    }

    // I.4: Compute network with strongly typed interface
    Graph graph;

    if (algorithm == "tcs") {
        graph = TCSAlgorithm::compute(haplotypes);
    } else if (algorithm == "msn") {
        graph = MSNAlgorithm::compute(haplotypes, input.epsilon);
    } else if (algorithm == "mjn") {
        graph = MJNAlgorithm::compute(haplotypes, input.epsilon);
    }

    // Extract output data
    for (const auto& vertex : graph.vertices()) {
        output.node_ids.push_back(vertex.id);
        output.node_sequences.push_back(vertex.sequence);
        output.node_samples.push_back(vertex.samples);
        output.node_is_median.push_back(vertex.is_median);
    }

    for (const auto& edge : graph.edges()) {
        output.edge_sources.push_back(edge.source);
        output.edge_targets.push_back(edge.target);
        output.edge_weights.push_back(edge.weight);
    }

    return output;
}

int NetworkBackend::run(int argc, char* argv[]) {
    try {
        // Parse algorithm
        std::string algorithm = parse_algorithm(argc, argv);
        if (algorithm.empty()) {
            return 1;
        }

        // Read input JSON
        std::string input_json = read_stdin();

        InputData input = parse_input_json(input_json);

        // Compute network
        OutputData output = compute_network(algorithm, input);

        // Write output JSON
        std::string output_json = generate_output_json(output);
        std::cout << output_json;

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}

} // namespace network
} // namespace haplokit
