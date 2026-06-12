// kcore_ldp_paper_style.h
#pragma once
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <iomanip>
#include <sstream>

#include "kcored/kcore_ldp.h"

// Graph结构体在graph.h中定义

// ===== Geometric sampler (port from reference coreness_estimation.cpp) =====
// Implements a two-sided geometric (discrete Laplace) sampler matching the reference
struct GeomDistribution {
    double lambda;
    std::mt19937_64 rng;
    std::uniform_real_distribution<double> uniform_dist;

    explicit GeomDistribution(double lambda_)
        : lambda(lambda_), rng(std::random_device{}()), uniform_dist(0.0, 1.0) {}

    int64_t geometric() {
        if (uniform_dist(rng) > -1.0 * std::expm1(-1.0 * lambda * std::numeric_limits<int64_t>::max())) {
            return std::numeric_limits<int64_t>::max();
        }

        int64_t left = 0;
        int64_t right = std::numeric_limits<int64_t>::max();

        while (left + 1 < right) {
            int64_t mid = left - static_cast<int64_t>(std::floor(
                (std::log(0.5) + std::log1p(std::exp(lambda * static_cast<double>(left - right)))) / lambda));

            if (mid <= left) {
                mid = left + 1;
            } else if (mid >= right) {
                mid = right - 1;
            }

            double q = std::expm1(lambda * static_cast<double>(left - mid)) /
                       std::expm1(lambda * static_cast<double>(left - right));

            if (uniform_dist(rng) <= q) {
                right = mid;
            } else {
                left = mid;
            }
        }

        return right;
    }

    int64_t twoSidedGeometric() {
        int64_t sample = 0;
        int64_t sign = -1;

        while (sample == 0 && sign == -1) {
            sample = geometric() - 1;
            sign = (uniform_dist(rng) < 0.5) ? -1 : 1;
        }

        return sample * sign;
    }
};

// ===== 参数与返回结构 =====
struct KCoreLDPParams {
    double epsilon = 1.0;      // 总隐私预算（对应YAML中的epsilon）
    double split = 0.8;        // 隐私预算拆分比例：ε1 = split*ε（阈值分组）, ε2 = (1-split)*ε（Level Moving）
    double phi = 0.5;          // 近似因子 φ（实验用，与非隐私k-core对齐）
    // 为与论文实验保持一致，设定 (2+η)=5.625 => η=3.625
    double eta = 3.625;          // 在 Go 实现中未直接使用，保留以便理论参数对齐
    bool   bias = true;        // 是否使用偏置校正（对应YAML中的bias）
    double bias_factor = 8.0;  // 偏置因子（对应YAML中的bias_factor）
    bool   noise = true;       // 是否添加噪声（对应YAML中的noise）
    int    num_workers = 80;   // 并行工作线程数（对应YAML中的num_workers）
    int    runs = 5;           // 重复运行次数（对应YAML中的runs）
    std::string partition_dir; // 分区图所在目录（例如 graphs/gplus_partitioned_80）
    bool distributed_mode = false; // true 表示启用 Go 风格的 worker/coord 模式
    unsigned int seed = 0;     // 可选随机种子
    bool aggregate_mean = true; // 多次运行是否取均值
    bool write_results = false; // 是否输出每次运行结果
    bool write_average = false; // 是否输出平均结果
    std::string result_dir;     // 结果输出目录
    std::string result_prefix;  // 输出文件前缀
    std::string result_tag;     // 可选标签
    int run_offset = 0;         // 运行编号偏移
    bool go_style_naming = false; // 采用Go论文的命名格式
    std::string graph_name;     // 用于Go风格命名
    std::string algo_name = "kcoreLDP";
    
    // 构造函数 - 匹配YAML配置和论文参数
        KCoreLDPParams(double eps = 1.0, double split_val = 0.8, double phi_val = 0.5, double eta_val = 3.625,
                                     bool bias_enabled = true, double bias_fac = 8.0, bool noise_enabled = true, 
                                     int workers = 80, int num_runs = 5) 
                : epsilon(eps), split(split_val), phi(phi_val), eta(eta_val), bias(bias_enabled), 
                    bias_factor(bias_fac), noise(noise_enabled), num_workers(workers), runs(num_runs) {}
};

// 工具：按论文设定 L=ceil(log n)/4，并定义 F(r)
// computeL now uses phi to match the reference: levels_per_group := ceil(log_{1+phi}(n)) / 4
static inline int computeL(int n, double phi) {
    // guard
    int nn = std::max(2, n);
    double log_base = std::log(static_cast<double>(nn)) / std::log(1.0 + phi);
    double val = std::ceil(log_base) / 4.0;
    int L = std::max(1, static_cast<int>(std::floor(val + 1e-12))); // floor to int, keep >=1
    return L;
}
static inline int F_of_round(int r, int L) {
    // 简化：把轮次 r 映到其所在的“组索引”（每组 L 个 level）
    return r / std::max(1, L);
}

// helper: log base conversion (used by EstimateCoreNumbers)
static inline double log_a_to_base_b(int a, double b) {
    return std::log2(static_cast<double>(a)) / std::log2(b);
}

// 为了与论文一致，最终 core 由协调端估计；留一个 stub 接口
// 真实实现应使用 LDS（各点最终层级/组）、L、lambda、psi 等来回推 core numbers
static std::vector<int> EstimateCoreNumbersStub(
    const std::vector<uint32_t>& final_level,
    int L, double lambda, double psi, double phi)
{
    // Ported coordinator mapping: compute levels_per_group then map levels -> core numbers
    const int n = static_cast<int>(final_level.size());
    std::vector<int> core(n, 0);

    // levels_per_group := ceil(log_a_to_base_b(n, 1.0+phi)) / 4 (use phi as in coordinator)
    const double levels_per_group = std::ceil(log_a_to_base_b(std::max(2, n), 1.0 + phi)) / 4.0;

    const double two_plus_lambda = 2.0 + lambda;
    const double one_plus_phi = 1.0 + phi;

    for (int i = 0; i < n; ++i) {
        double node_level = static_cast<double>(final_level[i]);
        double frac_numerator = node_level + 1.0;
        double power = std::max(std::floor(frac_numerator / levels_per_group) - 1.0, 0.0);
        double est = two_plus_lambda * std::pow(one_plus_phi, power);
        core[i] = static_cast<int>(std::round(est));
    }
    return core;
}

// ===== 主算法：论文 3.1 节结构化实现（分两阶段） =====
struct DegreeThresholdInfo {
    int threshold = 1;     // 参与轮次上限，同时作为 level-moving 阶段的分母
    bool permZero = false; // 永久跳过位：达到阈值或检查失败后标记
};

struct LDSLevelingState {
    std::vector<uint32_t> level;     // 每点当前层
    std::vector<int> group_index;    // F(r)：把轮次 r 映到组索引
    int num_rounds = 0;
    int L = 0;
};

// 阶段 1：Degree Thresholding（带噪度数、阈值 bias、得到每个点的轮次阈值）
// 返回每个点的 threshold，并通过引用参数给出全局最大阈值
static std::vector<DegreeThresholdInfo>
DegreeThresholding(const Graph& graph, const KCoreLDPParams& P, int /*L*/, int& global_round_threshold)
{
    const int n = graph.num_nodes();
    // 根据论文与参考实现：levels_per_group := ceil(log_a_to_base_b(n,1+phi))/4
    const double levels_per_group = std::ceil(log_a_to_base_b(std::max(2, n), 1.0 + P.phi)) / 4.0;

    // ε1 used as lambda in reference; use GeomDistribution sampler to match reference
    const double lambda1 = P.split * P.epsilon;
    GeomDistribution noise_deg(std::max(1e-12, lambda1/2.0));

    std::vector<DegreeThresholdInfo> info(n);
    int max_round_threshold = 0;

    for (int v=0; v<n; ++v) {
        const int d = int(graph.neighbor[v].size());
        int64_t noised_degree = int64_t(d);
        if (P.noise) {
            int64_t noise_sampled = noise_deg.twoSidedGeometric();
            noised_degree += noise_sampled;
            // bias subtraction as in reference: min(bias_factor * (2*exp(lambda))/(exp(2*lambda)-1), noised_degree)
            double cap = 0.0;
            if (lambda1 > 0) {
                cap = std::min( double(P.bias_factor) * (2.0 * std::exp(lambda1)) / (std::exp(2.0*lambda1) - 1.0), double(noised_degree) );
            }
            noised_degree = std::max<int64_t>(0, noised_degree - static_cast<int64_t>(std::floor(cap)) );
            // ensure degree at least 1 (matching reference adds 1)
            noised_degree = std::max<int64_t>(1, noised_degree + 1);
        }

        // threshold := ceil(log2(noised_degree)) * levels_per_group
        double thr = 1.0;
        if (noised_degree > 0) {
            thr = std::ceil(log_a_to_base_b(static_cast<int>(noised_degree), 2.0)) * levels_per_group;
        }
        // store as integer round-threshold (reference uses int(threshold)+1)
        info[v].threshold = static_cast<int>(thr) + 1;
        max_round_threshold = std::max(max_round_threshold, info[v].threshold);
    }
    global_round_threshold = max_round_threshold;
    return info;
}

// 阶段 2：Level Moving（同层邻居计数 + 几何噪声 + 计数侧 bias；与 (1+phi)^F(r) 比较）
static LDSLevelingState
LevelMoving(const Graph& graph, const KCoreLDPParams& P,
            std::vector<DegreeThresholdInfo>& info, int L, int max_round_threshold)
{
    const int n = graph.num_nodes();
    const double lambda2 = (1.0 - P.split) * P.epsilon;

    LDSLevelingState S;
    S.level.assign(n, 0);
    S.L = L;

    const double log_term = log_a_to_base_b(std::max(2, n), 1.0 + P.phi);
    const int analytic_rounds = static_cast<int>(std::ceil(4.0 * std::pow(log_term, 1.2)));
    const int cap_from_analysis = std::max(0, analytic_rounds - 2);
    const int round_limit = std::min(max_round_threshold, cap_from_analysis);
    S.num_rounds = std::max(0, round_limit);

    // group index per round: group = floor(round / levels_per_group)
    const double levels_per_group = std::ceil(log_a_to_base_b(std::max(2, n), 1.0 + P.phi)) / 4.0;
    S.group_index.resize(S.num_rounds);
    for (int r=0; r<S.num_rounds; ++r) S.group_index[r] = static_cast<int>(std::floor(double(r) / levels_per_group));

    for (int r=0; r<S.num_rounds; ++r) {
        const int Fr = S.group_index[r];
        std::vector<char> upgrade(n, 0);

        for (int v=0; v<n; ++v) {
            if (info[v].permZero) continue;
            if (r >= info[v].threshold) {
                info[v].permZero = true;
                continue;
            }
            if ((int)S.level[v] != r) continue;

            int U = 0;
            for (int u : graph.neighbor[v]) {
                if ((int)S.level[u] == r) ++U;
            }

            const int t = std::max(1, info[v].threshold);
            GeomDistribution g(std::max(1e-12, lambda2 / (2.0 * double(t))));

            int64_t U_tilde = int64_t(U);
            if (P.noise) {
                int64_t noise_sampled = g.twoSidedGeometric();
                U_tilde += noise_sampled;
                // extra_bias as in reference
                double scale = std::max(1e-12, lambda2 / (2.0 * double(t)));
                int64_t extra_bias = static_cast<int64_t>( std::floor(3.0 * (2.0 * std::exp(scale)) / std::pow((std::exp(2.0*scale) - 1.0), 3)) );
                U_tilde += extra_bias;
            }

            // upgrade condition matches Go implementation: compare against (1+phi)^{group}
            const double rhs = std::pow(1.0 + P.phi, static_cast<double>(std::max(0, Fr)));
            if (double(U_tilde) > rhs) {
                upgrade[v] = 1;
            } else {
                info[v].permZero = true;
            }
        }

        for (int v=0; v<n; ++v) if (upgrade[v]) S.level[v] += 1;
    }

    return S;
}

// ===== 外部主入口：与论文结构一致 =====
std::vector<int> estimateKCoreLDP_paper_style(const Graph& graph, const KCoreLDPParams& P,
                                              double lambda_for_stub = 0.0, double psi_for_stub = 0.0)
{
    if (P.distributed_mode && !P.partition_dir.empty()) {
        kcored::RunConfig config;
        config.n = graph.num_nodes();
        config.psi = P.phi;
        config.epsilon = P.epsilon;
        config.split = P.split;
        config.bias = P.bias;
        config.biasFactor = static_cast<int>(std::lround(P.bias_factor));
        config.noise = P.noise;
        config.numWorkers = P.num_workers;
        config.partitionDir = P.partition_dir;
        config.lambdaCoordinator = (lambda_for_stub != 0.0) ? lambda_for_stub : 0.5;
        config.writePerRun = P.write_results;
        config.writeAverage = P.write_average && P.aggregate_mean;
        config.resultDir = P.result_dir;
        config.resultPrefix = P.result_prefix;
        config.resultTag = P.result_tag;
        config.runOffset = P.run_offset;
        if (config.writePerRun && P.go_style_naming && !P.graph_name.empty()) {
            std::string graph = P.graph_name;
            std::string algo = P.algo_name.empty() ? std::string("kcoreLDP") : P.algo_name;
            int biasFlag = P.bias ? 1 : 0;
            int noiseFlag = P.noise ? 1 : 0;
            int biasFactor = static_cast<int>(std::lround(P.bias_factor));
            int workers = P.num_workers;
            double factor = P.split;
            double epsilon = P.epsilon;
            std::string tag = P.result_tag.empty() ? std::string("default") : P.result_tag;
            config.resultPrefix.clear();
            config.resultTag.clear();
            config.perRunFilenames.clear();
            config.perRunFilenames.reserve(std::max(1, P.runs));
            for (int run = 0; run < std::max(1, P.runs); ++run) {
                std::ostringstream oss;
                oss.setf(std::ios::fixed);
                oss << graph << "_" << algo << "_" << std::setprecision(2) << factor
                    << "_" << biasFlag << "_" << noiseFlag << "_" << biasFactor
                    << "_" << (P.run_offset + run) << "_" << workers << "_"
                    << std::setprecision(2) << epsilon << "_" << tag << ".txt";
                config.perRunFilenames.push_back(oss.str());
            }
            if (config.writeAverage) {
                std::ostringstream mean_ss;
                mean_ss.setf(std::ios::fixed);
                mean_ss << graph << "_" << algo << "_" << std::setprecision(2) << factor
                        << "_" << biasFlag << "_" << noiseFlag << "_" << biasFactor
                        << "_mean_" << workers << "_" << std::setprecision(2) << epsilon
                        << "_" << tag << ".txt";
                config.meanFilename = mean_ss.str();
            }
        }

        std::vector<double> estimates;
        if (P.aggregate_mean) {
            estimates = kcored::RunMultiple(config, std::max(1, P.runs), nullptr, P.seed);
        } else {
            auto single = kcored::RunSingle(config, P.seed);
            estimates = std::move(single.coreNumbers);
        }
        std::vector<int> rounded(estimates.size());
        for (size_t i = 0; i < estimates.size(); ++i) {
            rounded[i] = static_cast<int>(std::round(estimates[i]));
        }
        return rounded;
    }

    const int n = graph.num_nodes();
    const int L = computeL(n, P.phi);

    // 阶段 1：度阈值分组（带噪 & 阈值侧 bias）
    int max_round_threshold = 0;
    auto info = DegreeThresholding(graph, P, L, max_round_threshold);

    // 阶段 2：Level Moving（同层计数 + 几何噪声 + 计数侧 bias + 移动判据）
    auto S = LevelMoving(graph, P, info, L, max_round_threshold);

    // 最终 core 估计——交给协调端；这里放一个 stub
    // 论文：EstimateCoreNumbers(LDS, L, λ, ψ)
    return EstimateCoreNumbersStub(S.level, L, lambda_for_stub, psi_for_stub, P.phi);
}
