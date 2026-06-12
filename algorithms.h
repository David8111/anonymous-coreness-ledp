#pragma once

#include "graph.h"

#include <string>
#include <vector>

struct GlobalMetrics {
    double mae = 0.0;
    double rmse = 0.0;
    double acc_at_4 = 0.0;
    double mean_factor = 0.0;
    double p80_factor = 0.0;
    double p95_factor = 0.0;
};

struct AlgoResult {
    std::string name;
    std::vector<int> estimate;
    GlobalMetrics metrics;
    double runtime_seconds = 0.0;
};

std::vector<int> compute_coreness(Graph &g);
std::vector<int> estimate_coreness_LDP_hindex_iteration(Graph &g, double epsilon, int seed, int max_iter);
std::vector<int> estimate_coreness_LDP_hindex_transition(Graph &g, double epsilon, int seed, int max_iter);
std::vector<int> estimate_coreness_LDP_hindex_active(Graph &g, double epsilon, int seed, int max_iter);

bool compute_global_metrics(const std::vector<int> &true_vals,
                            const std::vector<int> &pred_vals,
                            int k,
                            GlobalMetrics &out);

std::vector<AlgoResult> run_four_algo_suite(Graph &g,
                                            double epsilon,
                                            const std::string &partition_dir,
                                            int total_rounds = 5,
                                            int seed = 42);

void print_global_summary_table(const std::vector<AlgoResult> &results);
bool write_metrics_csv(const std::vector<AlgoResult> &results, const std::string &path);
