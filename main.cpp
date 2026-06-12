#include "algorithms.h"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

int main(int argc, char *argv[])
{
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <epsilon> <graph_path> <partition_dir> [csv_output]\n";
        return 1;
    }

    const double epsilon = std::stod(argv[1]);
    const std::string graph_path = argv[2];
    const std::string partition_dir = argv[3];
    const std::string csv_output = (argc > 4) ? argv[4] : "";

    Graph g(graph_path, false);
    auto results = run_four_algo_suite(g, epsilon, partition_dir, 5, 42);

    std::cout << "nodes=" << g.num_nodes() << "\n";
    print_global_summary_table(results);

    if (!csv_output.empty()) {
        std::filesystem::create_directories(std::filesystem::path(csv_output).parent_path());
        if (!write_metrics_csv(results, csv_output)) {
            std::cerr << "Failed to write CSV: " << csv_output << "\n";
            return 2;
        }
        std::cout << "\nSaved metrics to " << csv_output << "\n";
    }

    return 0;
}
