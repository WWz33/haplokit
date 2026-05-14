#ifndef HAPLOKIT_NETWORK_BACKEND_H
#define HAPLOKIT_NETWORK_BACKEND_H

#include <string>
#include <vector>

namespace haplokit {
namespace network {

// JSON interface for network computation
// Reads JSON from stdin, writes JSON to stdout
class NetworkBackend {
public:
    // Main entry point
    static int run(int argc, char* argv[]);

private:
    // Parse command-line arguments
    static std::string parse_algorithm(int argc, char* argv[]);

    // Read JSON input from stdin
    static std::string read_stdin();

    // Parse input JSON
    struct InputData {
        std::vector<std::string> sequences;
        std::vector<std::vector<std::string>> samples;
        int epsilon = 0;
    };
    static InputData parse_input_json(const std::string& json);

    // Generate output JSON
    struct OutputData {
        std::vector<int> node_ids;
        std::vector<std::string> node_sequences;
        std::vector<std::vector<std::string>> node_samples;
        std::vector<bool> node_is_median;
        std::vector<int> edge_sources;
        std::vector<int> edge_targets;
        std::vector<int> edge_weights;
    };
    static std::string generate_output_json(const OutputData& data);

    // Compute network using specified algorithm
    static OutputData compute_network(const std::string& algorithm,
                                      const InputData& input);
};

} // namespace network
} // namespace haplokit

#endif // HAPLOKIT_NETWORK_BACKEND_H
