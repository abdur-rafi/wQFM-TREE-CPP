// Micro-benchmark of the scorer kernels: score() and gains() on a top-level subproblem (all taxa real)
// and on subproblems with dummy taxa, with a fixed random bipartition. Prints the time per call (wall
// and CPU: the CPU time is not affected by the machine sleeping), ns per forest node, and a checksum of
// all values, so kernel changes can be checked for bit-identical results.
//
// usage: bench_scorer <gene trees> [threads=1] [repetitions=20] [weighting=none|support|length|hybrid]
//   weighting other than none benchmarks the weighted scorer (WScorer; supports on a 0..1 scale)

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "forest.hpp"
#include "newick.hpp"
#include "scorer.hpp"
#include "thread_pool.hpp"
#include "units.hpp"
#include "wscorer.hpp"

#include <sys/resource.h>

using namespace wqfm;

namespace {

std::vector<ParsedTree> read_trees(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read " + path);
    std::vector<ParsedTree> trees;
    std::string line;
    while (std::getline(in, line))
        if (line.find('(') != std::string::npos) trees.push_back(parse_newick(line));
    return trees;
}

uint64_t bits(double x) {
    uint64_t u;
    std::memcpy(&u, &x, sizeof u);
    return u;
}

double cpu_seconds() {
    rusage u;
    getrusage(RUSAGE_SELF, &u);
    return u.ru_utime.tv_sec + u.ru_stime.tv_sec + 1e-6 * (u.ru_utime.tv_usec + u.ru_stime.tv_usec);
}

void bench(const char* name, IScorer& sc, const Subproblem& sp, int reps) {
    sc.set_partition(sp.part_r, sp.part_d);
    std::vector<double> gr, gd;
    uint64_t h = 0;
    auto mix = [&](double x) { h = (h ^ bits(x)) * 0x100000001B3ULL; };

    using Clock = std::chrono::steady_clock;
    auto t0 = Clock::now();
    double c0 = cpu_seconds();
    for (int r = 0; r < reps; ++r) mix(sc.score());
    const double ts = std::chrono::duration<double>(Clock::now() - t0).count() / reps;
    const double cs = (cpu_seconds() - c0) / reps;
    t0 = Clock::now();
    c0 = cpu_seconds();
    for (int r = 0; r < reps; ++r) {
        mix(sc.gains(gr, gd));
        for (double x : gr) mix(x);
        for (double x : gd) mix(x);
    }
    const double tg = std::chrono::duration<double>(Clock::now() - t0).count() / reps;
    const double cg = (cpu_seconds() - c0) / reps;
    const int N = sc.forest().nodes();
    std::printf("%-8s nodes %8d  score %8.3f ms (cpu %8.3f, %6.2f ns/node)  gains %9.3f ms (cpu %9.3f, %7.2f ns/node)  checksum %016llx\n",
                name, N, ts * 1e3, cs * 1e3, cs * 1e9 / N, tg * 1e3, cg * 1e3, cg * 1e9 / N, static_cast<unsigned long long>(h));
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: bench_scorer <gene trees> [threads] [repetitions] [none|support|length|hybrid]\n");
        return 2;
    }
    const int threads = argc > 2 ? std::atoi(argv[2]) : 1;
    const int reps = argc > 3 ? std::atoi(argv[3]) : 20;
    const std::string weighting = argc > 4 ? argv[4] : "none";
    if (weighting != "none" && weighting != "support" && weighting != "length" && weighting != "hybrid") {
        std::fprintf(stderr, "unknown weighting %s\n", weighting.c_str());
        return 2;
    }
    const std::vector<ParsedTree> trees = read_trees(argv[1]);
    Taxa taxa;
    for (const auto& t : trees)
        for (size_t v = 0; v < t.label.size(); ++v)
            if (t.children[v].empty()) taxa.add(t.label[v]);
    const Forest forest = make_forest(trees, taxa);
    const int n = forest.n;
    std::unique_ptr<ThreadPool> pool;
    if (threads > 1) pool = std::make_unique<ThreadPool>(threads);
    std::mt19937_64 rng(42);
    auto make = [&]() -> std::unique_ptr<IScorer> {
        if (weighting == "none") return std::make_unique<Scorer>(forest);
        WScorer::Options o;
        o.support = weighting == "support" || weighting == "hybrid";
        o.length = weighting == "length" || weighting == "hybrid";
        return std::make_unique<WScorer>(forest, true, o);
    };

    // top level: every taxon real
    Subproblem top;
    DummyPool dp;
    dp.reserve(64);
    for (int t = 0; t < n; ++t) {
        top.reals.push_back(t);
        top.part_r.push_back(static_cast<int8_t>(rng() & 1));
    }
    std::unique_ptr<IScorer> s1p = make();
    IScorer& s1 = *s1p;
    s1.set_pool(pool.get());
    s1.setup(top, dp);
    bench("top", s1, top, reps);

    // half the taxa real, the other half one dummy taxon (a second-level subproblem)
    Subproblem half;
    std::vector<int> other;
    for (int t = 0; t < n; ++t) {
        if (t % 2 == 0) {
            half.reals.push_back(t);
            half.part_r.push_back(static_cast<int8_t>(rng() & 1));
        } else {
            other.push_back(t);
        }
    }
    half.dummies.push_back(dp.add(other, {}));
    half.part_d.push_back(0);
    std::unique_ptr<IScorer> s3p = make();
    IScorer& s3 = *s3p;
    s3.set_pool(pool.get());
    s3.setup(half, dp);
    bench("1 dummy", s3, half, reps);

    // a third of the taxa real, the rest under 4 dummy taxa
    Subproblem sub;
    std::vector<std::vector<int>> groups(4);
    for (int t = 0; t < n; ++t) {
        if (t % 3 == 0) {
            sub.reals.push_back(t);
            sub.part_r.push_back(static_cast<int8_t>(rng() & 1));
        } else {
            groups[rng() % 4].push_back(t);
        }
    }
    for (auto& gset : groups) {
        sub.dummies.push_back(dp.add(gset, {}));
        sub.part_d.push_back(static_cast<int8_t>(rng() & 1));
    }
    std::unique_ptr<IScorer> s2p = make();
    IScorer& s2 = *s2p;
    s2.set_pool(pool.get());
    s2.setup(sub, dp);
    bench("4 dummy", s2, sub, reps);

    // deep subproblem: a tenth of the taxa real, 30 small dummy taxa and one for the rest of the tree
    Subproblem deep;
    std::vector<int> order(n);
    for (int t = 0; t < n; ++t) order[t] = t;
    std::shuffle(order.begin(), order.end(), rng);
    const int nr = n / 10, small = std::max(1, n / 40);
    int at = 0;
    for (; at < nr; ++at) {
        deep.reals.push_back(order[at]);
        deep.part_r.push_back(static_cast<int8_t>(rng() & 1));
    }
    for (int k = 0; k < 30 && at + small < n; ++k, at += small) {
        deep.dummies.push_back(dp.add(std::vector<int>(order.begin() + at, order.begin() + at + small), {}));
        deep.part_d.push_back(static_cast<int8_t>(rng() & 1));
    }
    deep.dummies.push_back(dp.add(std::vector<int>(order.begin() + at, order.end()), {}));
    deep.part_d.push_back(static_cast<int8_t>(rng() & 1));
    std::unique_ptr<IScorer> s4p = make();
    IScorer& s4 = *s4p;
    s4.set_pool(pool.get());
    s4.setup(deep, dp);
    bench("31 dummy", s4, deep, reps);
    return 0;
}
