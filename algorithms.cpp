#include "algorithms.h"

#include "single_coreness_function.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <future>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <random>

using namespace std;

template <typename Func>
static void parallel_for(int start, int end, Func func, int num_threads = std::thread::hardware_concurrency())
{
    if (num_threads <= 1 || end <= start) {
        for (int i = start; i < end; ++i) func(i);
        return;
    }

    std::vector<std::future<void>> futures;
    int chunk_size = std::max(1, (end - start) / num_threads);
    for (int t = 0; t < num_threads; ++t) {
        int thread_start = start + t * chunk_size;
        if (thread_start >= end) break;
        int thread_end = (t == num_threads - 1) ? end : std::min(end, start + (t + 1) * chunk_size);
        futures.push_back(std::async(std::launch::async, [=, &func]() {
            for (int i = thread_start; i < thread_end; ++i) func(i);
        }));
    }
    for (auto &future : futures) future.wait();
}

static int sample_symmetric_geometric(double epsilon, std::mt19937 &rng)
{
    if (!(epsilon > 0.0) || !std::isfinite(epsilon)) return 0;
    std::geometric_distribution<int> geo(1.0 - std::exp(-epsilon));
    return geo(rng) - geo(rng);
}

static double discrete_laplace_expected_abs(double epsilon)
{
    if (!(epsilon > 0.0) || !std::isfinite(epsilon)) {
        return std::numeric_limits<double>::infinity();
    }
    const double a = std::exp(-epsilon);
    return (2.0 * a) / std::max(1e-12, 1.0 - a * a);
}

static int compute_h_index(const std::vector<int> &degs)
{
    if (degs.empty()) return 0;
    std::vector<int> sorted = degs;
    std::sort(sorted.begin(), sorted.end(), std::greater<int>());

    int h = 0;
    for (int d : sorted) {
        if (d >= h + 1) ++h;
        else break;
    }
    return h;
}

static std::vector<int> initialize_noisy_degree_round(Graph &g, double epsilon_round, std::mt19937 &rng)
{
    const int n = g.num_nodes();
    std::vector<int> noisy_degree(n, 0);
    for (int u = 0; u < n; ++u) {
        int deg = static_cast<int>(g.neighbor[u].size());
        noisy_degree[u] = std::max(0, deg + sample_symmetric_geometric(epsilon_round, rng));
    }
    return noisy_degree;
}

static std::vector<int> release_noisy_hindex_round(Graph &g,
                                                   const std::vector<int> &state,
                                                   double epsilon_round,
                                                   std::mt19937 &rng,
                                                   const std::vector<char> *active = nullptr)
{
    const int n = g.num_nodes();
    std::vector<int> next = state;
    if (!(epsilon_round > 0.0)) return next;

    std::vector<int> noise_values(n, 0);
    for (int u = 0; u < n; ++u) noise_values[u] = sample_symmetric_geometric(epsilon_round, rng);

    parallel_for(0, n, [&](int u) {
        if (active && !(*active)[u]) return;
        std::vector<int> neigh_state;
        neigh_state.reserve(g.neighbor[u].size());
        for (int v : g.neighbor[u]) neigh_state.push_back(state[v]);
        next[u] = std::max(0, compute_h_index(neigh_state) + noise_values[u]);
    });
    return next;
}

static double mean_abs_delta(const std::vector<int> &prev,
                             const std::vector<int> &next,
                             const std::vector<char> *mask = nullptr)
{
    if (prev.size() != next.size() || prev.empty()) return 0.0;
    long double sum = 0.0L;
    int cnt = 0;
    for (size_t i = 0; i < prev.size(); ++i) {
        if (mask && !(*mask)[i]) continue;
        sum += std::abs(next[i] - prev[i]);
        ++cnt;
    }
    return cnt ? static_cast<double>(sum / static_cast<long double>(cnt)) : 0.0;
}

static int count_active_nodes(const std::vector<char> &active)
{
    return static_cast<int>(std::count(active.begin(), active.end(), static_cast<char>(1)));
}

static double state_entropy_bits(const std::vector<int> &state)
{
    if (state.empty()) return 0.0;
    std::vector<int> sorted = state;
    std::sort(sorted.begin(), sorted.end());
    const double n = static_cast<double>(sorted.size());
    double entropy = 0.0;
    int run = 1;
    for (size_t i = 1; i <= sorted.size(); ++i) {
        if (i < sorted.size() && sorted[i] == sorted[i - 1]) {
            ++run;
            continue;
        }
        const double p = static_cast<double>(run) / n;
        entropy -= p * std::log2(std::max(1e-12, p));
        run = 1;
    }
    return entropy;
}

static double state_entropy_bits_masked(const std::vector<int> &state, const std::vector<char> &mask)
{
    if (state.empty() || state.size() != mask.size()) return 0.0;
    std::vector<int> vals;
    vals.reserve(state.size());
    for (size_t i = 0; i < state.size(); ++i) {
        if (mask[i]) vals.push_back(state[i]);
    }
    return vals.empty() ? 0.0 : state_entropy_bits(vals);
}

static std::vector<int> final_hindex_polish(Graph &g,
                                            const std::vector<int> &state,
                                            double epsilon_rem,
                                            std::mt19937 &rng)
{
    if (!(epsilon_rem > 0.0)) return state;
    return release_noisy_hindex_round(g, state, epsilon_rem, rng, nullptr);
}

std::vector<int> estimate_coreness_LDP_hindex_iteration(Graph &g, double epsilon, int seed, int max_iter)
{
    const int n = g.num_nodes();
    std::mt19937 rng(seed);

    const int total_rounds = std::max(1, max_iter);
    const int iter_rounds = total_rounds - 1;
    const double eps_round = epsilon / static_cast<double>(total_rounds);

    std::vector<int> coreness = initialize_noisy_degree_round(g, eps_round, rng);
    if (iter_rounds <= 0) return coreness;

    std::vector<int> noise_values(n * iter_rounds);
    for (int &value : noise_values) value = sample_symmetric_geometric(eps_round, rng);

    for (int it = 0; it < iter_rounds; ++it) {
        std::vector<int> next(n, 0);
        parallel_for(0, n, [&](int u) {
            std::vector<int> neigh_cor;
            neigh_cor.reserve(g.neighbor[u].size());
            for (int v : g.neighbor[u]) neigh_cor.push_back(coreness[v]);
            next[u] = std::max(0, compute_h_index(neigh_cor) + noise_values[it * n + u]);
        });
        coreness.swap(next);
    }
    return coreness;
}

std::vector<int> estimate_coreness_LDP_hindex_transition(Graph &g, double epsilon, int seed, int max_iter)
{
    std::mt19937 rng(seed);
    if (!(epsilon > 0.0)) return initialize_noisy_degree_round(g, 0.0, rng);

    const int total_rounds = std::max(1, max_iter);
    const int iter_rounds = total_rounds - 1;
    const double eps_round = epsilon / static_cast<double>(total_rounds);

    std::vector<int> state = initialize_noisy_degree_round(g, eps_round, rng);
    if (iter_rounds <= 0) return state;

    std::vector<char> active(g.num_nodes(), 1);
    double prev_eps = eps_round;

    for (int round = 0; round < iter_rounds; ++round) {
        if (count_active_nodes(active) == 0) break;
        std::vector<int> next = release_noisy_hindex_round(g, state, eps_round, rng, &active);
        const double noise_floor =
            discrete_laplace_expected_abs(prev_eps) + discrete_laplace_expected_abs(eps_round);

        for (int u = 0; u < g.num_nodes(); ++u) {
            if (active[u] && std::abs(next[u] - state[u]) <= noise_floor) active[u] = 0;
        }

        state.swap(next);
        prev_eps = eps_round;
    }

    return state;
}

std::vector<int> estimate_coreness_LDP_hindex_active(Graph &g, double epsilon, int seed, int max_iter)
{
    (void)max_iter;
    std::mt19937 rng(seed);
    if (!(epsilon > 0.0)) return initialize_noisy_degree_round(g, 0.0, rng);

    const int safety_rounds = std::max(1, g.num_nodes());
    const double init_entropy =
        std::max(1.0, std::log2(1.0 + static_cast<double>(std::max(2, g.num_nodes()))));
    const double eps_init = epsilon / (1.0 + init_entropy);

    std::vector<int> state = initialize_noisy_degree_round(g, eps_init, rng);
    std::vector<char> active(g.num_nodes(), 1);
    double spent = eps_init;
    double prev_eps = eps_init;

    for (int round = 1; round < safety_rounds; ++round) {
        if (count_active_nodes(active) == 0) break;
        const double remaining = std::max(0.0, epsilon - spent);
        if (!(remaining > 0.0)) break;

        const double entropy = std::max(0.0, state_entropy_bits_masked(state, active));
        const double eps_t = remaining / (1.0 + entropy);
        if (!(eps_t > 1e-6)) break;

        std::vector<int> next = release_noisy_hindex_round(g, state, eps_t, rng, &active);
        spent += eps_t;
        const double delta = mean_abs_delta(state, next, &active);
        const double noise_floor =
            discrete_laplace_expected_abs(prev_eps) + discrete_laplace_expected_abs(eps_t);

        for (int u = 0; u < g.num_nodes(); ++u) {
            if (active[u] && std::abs(next[u] - state[u]) <= noise_floor) active[u] = 0;
        }

        state.swap(next);
        prev_eps = eps_t;
        if (delta <= noise_floor || count_active_nodes(active) == 0) break;
    }

    return final_hindex_polish(g, state, std::max(0.0, epsilon - spent), rng);
}

std::vector<int> compute_coreness(Graph &g)
{
    const int n = g.num_nodes();
    std::vector<int> degree(n);
    for (int v = 0; v < n; ++v) degree[v] = static_cast<int>(g.neighbor[v].size());

    int max_deg = *std::max_element(degree.begin(), degree.end());
    std::vector<int> bin(max_deg + 1, 0);
    for (int d : degree) ++bin[d];

    int start = 0;
    for (int d = 0; d <= max_deg; ++d) {
        int cnt = bin[d];
        bin[d] = start;
        start += cnt;
    }

    std::vector<int> vert(n), pos(n);
    for (int v = 0; v < n; ++v) {
        int d = degree[v];
        vert[bin[d]] = v;
        pos[v] = bin[d];
        bin[d]++;
    }

    for (int d = max_deg; d > 0; --d) bin[d] = bin[d - 1];
    bin[0] = 0;

    std::vector<int> coreness(n);
    for (int i = 0; i < n; ++i) {
        int v = vert[i];
        coreness[v] = degree[v];
        for (int u : g.neighbor[v]) {
            if (degree[u] > degree[v]) {
                int du = degree[u];
                int pu = pos[u];
                int pw = bin[du];
                int w = vert[pw];
                if (u != w) {
                    vert[pu] = w;
                    pos[w] = pu;
                    vert[pw] = u;
                    pos[u] = pw;
                }
                bin[du]++;
                degree[u]--;
            }
        }
    }
    return coreness;
}

bool compute_global_metrics(const std::vector<int> &true_vals,
                            const std::vector<int> &pred_vals,
                            int k,
                            GlobalMetrics &out)
{
    if (true_vals.size() != pred_vals.size() || true_vals.empty()) return false;

    const size_t n = true_vals.size();
    double mae = 0.0;
    double mse = 0.0;
    size_t agree_k = 0;
    std::vector<double> factors;
    factors.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        const double err = std::abs(true_vals[i] - pred_vals[i]);
        mae += err;
        mse += err * err;
        if ((true_vals[i] >= k) == (pred_vals[i] >= k)) ++agree_k;

        int true_c = std::max(true_vals[i], 1);
        int pred_c = std::max(pred_vals[i], 1);
        factors.push_back(static_cast<double>(std::max(pred_c, true_c)) /
                          static_cast<double>(std::min(pred_c, true_c)));
    }

    mae /= static_cast<double>(n);
    mse /= static_cast<double>(n);
    std::sort(factors.begin(), factors.end());

    double sum = 0.0;
    size_t cnt = 0;
    for (double value : factors) {
        if (!std::isinf(value)) {
            sum += value;
            ++cnt;
        }
    }

    auto percentile = [&](double q) {
        size_t idx = static_cast<size_t>(std::ceil(q * factors.size()));
        if (idx) --idx;
        return factors[idx];
    };

    out.mae = mae;
    out.rmse = std::sqrt(mse);
    out.acc_at_4 = static_cast<double>(agree_k) / static_cast<double>(n);
    out.mean_factor = cnt ? sum / static_cast<double>(cnt) : 0.0;
    out.p80_factor = percentile(0.80);
    out.p95_factor = percentile(0.95);
    return true;
}

std::vector<AlgoResult> run_four_algo_suite(Graph &g,
                                            double epsilon,
                                            const std::string &partition_dir,
                                            int total_rounds,
                                            int seed)
{
    const std::vector<int> true_coreness = compute_coreness(g);
    std::vector<AlgoResult> results;
    results.reserve(4);

    auto run_one = [&](const std::string &name, auto fn) {
        auto t0 = std::chrono::high_resolution_clock::now();
        std::vector<int> pred = fn();
        auto t1 = std::chrono::high_resolution_clock::now();
        AlgoResult result;
        result.name = name;
        result.estimate = std::move(pred);
        result.runtime_seconds =
            std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1e6;
        compute_global_metrics(true_coreness, result.estimate, 4, result.metrics);
        results.push_back(std::move(result));
    };

    run_one("kcoreH", [&]() {
        return estimate_coreness_LDP_hindex_iteration(g, epsilon, seed, total_rounds);
    });

    run_one("kcored", [&]() {
        KCoreLDPParams params;
        params.epsilon = epsilon;
        params.split = 0.8;
        params.phi = 0.5;
        params.eta = 3.625;
        params.bias = true;
        params.bias_factor = 8.0;
        params.noise = true;
        params.num_workers = 80;
        params.partition_dir = partition_dir;
        params.distributed_mode = !partition_dir.empty();
        params.seed = static_cast<unsigned int>(seed);
        params.runs = 1;
        params.aggregate_mean = false;
        return estimateKCoreLDP_paper_style(g, params);
    });

    run_one("kcoreHT", [&]() {
        return estimate_coreness_LDP_hindex_transition(g, epsilon, seed, total_rounds);
    });

    run_one("kcoreHA", [&]() {
        return estimate_coreness_LDP_hindex_active(g, epsilon, seed, 0);
    });

    return results;
}

void print_global_summary_table(const std::vector<AlgoResult> &results)
{
    std::cout << "\n=== Global Metrics on Facebook ===\n\n";
    std::cout << std::left
              << std::setw(12) << "Algorithm"
              << std::setw(10) << "MAE"
              << std::setw(10) << "RMSE"
              << std::setw(10) << "Acc@4"
              << std::setw(14) << "MeanFactor"
              << std::setw(10) << "P80"
              << std::setw(10) << "P95"
              << std::setw(10) << "Time(s)" << "\n";
    std::cout << std::string(86, '-') << "\n";

    for (const auto &result : results) {
        std::cout << std::left
                  << std::setw(12) << result.name
                  << std::setw(10) << std::fixed << std::setprecision(4) << result.metrics.mae
                  << std::setw(10) << std::fixed << std::setprecision(4) << result.metrics.rmse
                  << std::setw(10) << std::fixed << std::setprecision(2) << (100.0 * result.metrics.acc_at_4)
                  << std::setw(14) << std::fixed << std::setprecision(4) << result.metrics.mean_factor
                  << std::setw(10) << std::fixed << std::setprecision(4) << result.metrics.p80_factor
                  << std::setw(10) << std::fixed << std::setprecision(4) << result.metrics.p95_factor
                  << std::setw(10) << std::fixed << std::setprecision(4) << result.runtime_seconds
                  << "\n";
    }
}

bool write_metrics_csv(const std::vector<AlgoResult> &results, const std::string &path)
{
    std::ofstream fout(path);
    if (!fout.is_open()) return false;

    fout << "algo,mae,rmse,acc_at_4,mean_factor,p80_factor,p95_factor,time_s\n";
    fout << std::fixed << std::setprecision(8);
    for (const auto &result : results) {
        fout << result.name << ","
             << result.metrics.mae << ","
             << result.metrics.rmse << ","
             << result.metrics.acc_at_4 << ","
             << result.metrics.mean_factor << ","
             << result.metrics.p80_factor << ","
             << result.metrics.p95_factor << ","
             << result.runtime_seconds << "\n";
    }
    return true;
}
