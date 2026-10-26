// wqfm-tree: species tree estimation from gene trees with the wQFM-TREE algorithm.
//
// usage: wqfm-tree -i <gene trees> -o <output> [-c <consensus tree>] [-t threads] [--relative-eps]
//                  [--resolve-polytomies] [--taxa-order <file>] [--trace <file>] [--replay <java trace>]
//                  [--subproblem-log <file>]
//
//   -i            gene trees, one Newick tree per line (branch lengths and support values are ignored)
//   -c            greedy consensus of the gene trees (e.g. PAUP's, as in run.sh); default: built in
//                 (complete gene trees only; PAUP may break frequency ties differently)
//   --write-consensus  write the consensus tree used
//   --resolve-polytomies  resolve gene-tree polytomies arbitrarily before the search, as
//                 runResolvingPolytomy.sh did (DendroPy: join the first two children until two remain).
//                 The built-in consensus still uses the unresolved trees, as that script's PAUP run did.
//   --taxa-order  "id label" lines, or a Java trace ("T id label" lines). Taxon ids decide ties, and
//                 the Java version numbers taxa in HashSet order; pass its order to reproduce it.
//                 Default: order of first appearance in the gene trees.
//   --trace       write a step-by-step trace (format shared with the instrumented Java build)
//   --replay      follow the decisions recorded in a Java trace (validation; implies its taxa order)
//   --subproblem-log  one tab-separated line per subproblem: real and dummy taxa, score of the initial
//                 bipartition, score after FM, FM passes, 1 if FM changed the bipartition (analysis)
//   -t            threads (default 1); the tree does not depend on it. Every FM move ends with a
//                 parallel loop, so keep the threads of concurrent jobs within the core count.
//   --quiet        no progress messages (stderr); --progress N prints the search progress every N seconds
//   --relative-eps FM keeps a longer prefix only if it beats the best by 1e-5 + 1e-12 * (sum of
//                 |gain| so far) instead of Java's 1e-5 (robust to rounding at large score scales;
//                 may differ from the Java version)
//   --fm-patience K  end an FM pass after K moves in a row without a new best prefix instead of moving
//                 every unit (faster; may differ from the Java version; 0: off, the default)
//   --fm-keep-flips  count an FM pass whose best prefix moves every unit as an improvement and start
//                 another pass, as Java does. By default it ends the subproblem's FM: it is the same
//                 bipartition with the sides swapped, its gain is rounding (same trees, faster)
//   --no-compress score on the full gene trees instead of the per-subproblem compressed ones (testing)
//   --kappa K, --power P, --dummy-gamma G  experimental objectives: maximise (2F - K T) / T^P, where
//                 F = satisfied + unresolved/2 and T = all quartets with two units per side (K = 1, P = 0:
//                 satisfied - violated; K = 2a/(1+a): satisfied - a violated), and weight the real taxa
//                 under a dummy taxon of m members by m^(G-1) per level (G = 0: total weight 1);
//                 --init-kappa K scores the initial-bipartition candidates with K instead
//   --weight support|length|hybrid  weighted quartets (as weighted ASTRAL): by 1 - prod(1 - s) over the
//                 supports s on the path between a quartet's anchors, by exp(-total pendant length), or
//                 both; --support-scale MIN MAX rescales supports from [MIN, MAX] (default: 0..1, or
//                 0..100 if any support exceeds 1)

#include <cctype>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "consensus.hpp"
#include "forest.hpp"
#include "greedy.hpp"
#include "newick.hpp"
#include "trace.hpp"
#include "wqfm.hpp"

using namespace wqfm;

namespace {

std::vector<ParsedTree> read_trees(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read " + path);
    std::vector<ParsedTree> trees;
    std::string line;
    for (int no = 1; std::getline(in, line); ++no) {
        if (line.find('(') == std::string::npos) continue;
        try {
            trees.push_back(parse_newick(line));
        } catch (const std::exception& e) {
            throw std::runtime_error(path + " line " + std::to_string(no) + ": " + e.what());
        }
    }
    if (trees.empty()) throw std::runtime_error("no trees in " + path);
    return trees;
}

Taxa taxa_from_order(const std::vector<std::pair<int, std::string>>& entries) {
    std::vector<std::string> by_id(entries.size());
    for (const auto& [id, label] : entries) {
        if (id < 0 || id >= static_cast<int>(entries.size()) || !by_id[id].empty())
            throw std::runtime_error("taxa order: ids must be 0..n-1, each once");
        by_id[id] = label;
    }
    Taxa taxa;
    for (const auto& l : by_id) taxa.add(l);
    return taxa;
}

std::vector<std::pair<int, std::string>> read_order(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read " + path);
    std::vector<std::pair<int, std::string>> out;
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream ss(line);
        std::string first;
        if (!(ss >> first)) continue;
        if (first == "T") { if (!(ss >> first)) continue; }
        else if (line.size() > 1 && line[1] == ' ' && std::isalpha(static_cast<unsigned char>(line[0]))) continue;   // other trace lines
        std::string label;
        if (!(ss >> label)) continue;
        out.emplace_back(std::stoi(first), label);
    }
    return out;
}

[[noreturn]] void usage() {
    std::fprintf(stderr, "usage: wqfm-tree -i <gene trees> -o <output> [-c <consensus tree>] [--write-consensus <file>]\n"
                         "                 [-t threads (default 1)] [--relative-eps] [--fm-patience K] [--fm-keep-flips] [--resolve-polytomies] [--taxa-order <file>] [--trace <file>]\n"
                         "                 [--replay <java trace>] [--no-compress] [--quiet] [--progress <seconds, default 10>]\n"
                         "                 [--subproblem-log <file>] [--weight support|length|hybrid] [--support-scale MIN MAX]\n");
    std::exit(2);
}

}  // namespace

int main(int argc, char** argv) {
    std::string in_path, cons_path, out_path, order_path, trace_path, replay_path, cons_out_path, sublog_path;
    SearchOptions opt;
    bool quiet = false, resolve = false;
    double support_max = -1;                     // --support-scale MAX; < 0: from the data
    opt.progress_secs = 10;
    for (int a = 1; a < argc; ++a) {
        const std::string arg = argv[a];
        auto value = [&]() -> std::string { if (a + 1 >= argc) usage(); return argv[++a]; };
        if (arg == "-i") in_path = value();
        else if (arg == "-c") cons_path = value();
        else if (arg == "-o") out_path = value();
        else if (arg == "--taxa-order") order_path = value();
        else if (arg == "--trace") trace_path = value();
        else if (arg == "--replay") replay_path = value();
        else if (arg == "--no-compress") opt.compress = false;
        else if (arg == "--relative-eps") opt.relative_eps = true;
        else if (arg == "--fm-patience") opt.fm_patience = std::max(0, std::atoi(value().c_str()));
        else if (arg == "--fm-keep-flips") opt.fm_skip_flips = false;
        else if (arg == "--resolve-polytomies") resolve = true;
        else if (arg == "--quiet") quiet = true;
        else if (arg == "--progress") opt.progress_secs = std::max(0.1, std::atof(value().c_str()));
        else if (arg == "--write-consensus") cons_out_path = value();
        else if (arg == "--subproblem-log") sublog_path = value();
        else if (arg == "--kappa") opt.kappa = std::atof(value().c_str());
        else if (arg == "--power") opt.power = std::atof(value().c_str());
        else if (arg == "--dummy-gamma") opt.dummy_gamma = std::atof(value().c_str());
        else if (arg == "--init-kappa") opt.init_kappa = std::atof(value().c_str());
        else if (arg == "--weight") {
            const std::string w = value();
            opt.weighting = w == "support" ? 1 : w == "length" ? 2 : w == "hybrid" ? 3 : -1;
            if (opt.weighting < 0) usage();
        }
        else if (arg == "--support-scale") { opt.support_min = std::atof(value().c_str()); support_max = std::atof(value().c_str()); }
        else if (arg == "-t") opt.threads = std::max(1, std::atoi(value().c_str()));
        else usage();
    }
    if (in_path.empty() || out_path.empty()) usage();
    if (quiet) opt.progress_secs = 0;

    try {
        const auto t0 = std::chrono::steady_clock::now();
        auto log = [&](const std::string& msg) {
            if (quiet) return;
            std::fprintf(stderr, "[%7.1f s] %s\n", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(), msg.c_str());
        };
        log("reading gene trees from " + in_path);
        std::vector<ParsedTree> gene_trees = read_trees(in_path);

        std::unique_ptr<Replay> replay;
        if (!replay_path.empty()) replay = std::make_unique<Replay>(Replay::load(replay_path));

        Taxa taxa;
        if (!order_path.empty()) taxa = taxa_from_order(read_order(order_path));
        else if (replay) taxa = taxa_from_order(replay->taxa);
        else for (const auto& t : gene_trees) for (size_t v = 0; v < t.label.size(); ++v) if (t.children[v].empty()) taxa.add(t.label[v]);

        Forest forest = make_forest(gene_trees, taxa);
        std::vector<ParsedTree>().swap(gene_trees);   // only the forest is used from here on (~60 bytes a node)
        log(std::to_string(forest.k) + " gene trees, " + std::to_string(taxa.size()) + " taxa, " + std::to_string(forest.nodes()) + " nodes");
        if (opt.weighting & 1) {
            double smax = 0;
            long with = 0, internal = 0;
            for (int v = 0; v < forest.nodes(); ++v) {
                if (forest.parent[v] < 0 || forest.child_off[v + 1] == forest.child_off[v]) continue;
                ++internal;
                if (!std::isnan(forest.supp[v])) { ++with; smax = std::max(smax, forest.supp[v]); }
            }
            if (support_max < 0) support_max = smax > 1 ? 100 : 1;
            opt.support_max = support_max;
            log(std::to_string(with) + " of " + std::to_string(internal) + " internal branches have supports (max " +
                std::to_string(smax) + "); rescaling from [" + std::to_string(opt.support_min) + ", " + std::to_string(opt.support_max) + "]");
        }
        std::vector<char> seen(taxa.size(), 0);
        for (int t : forest.taxon) if (t >= 0) seen[t] = 1;
        for (int t = 0; t < taxa.size(); ++t)
            if (!seen[t]) throw std::runtime_error("taxon '" + taxa.label[t] + "' from the taxa order is in no gene tree");

        std::string cons_newick;
        if (!cons_path.empty()) {
            std::ifstream in(cons_path);
            std::string line;
            while (std::getline(in, line) && line.find('(') == std::string::npos) {}
            if (line.find('(') == std::string::npos) throw std::runtime_error("no tree in " + cons_path);
            cons_newick = line;
        } else {
            log("building the greedy consensus tree");
            cons_newick = greedy_consensus_newick(forest, taxa);
        }
        if (!cons_out_path.empty()) std::ofstream(cons_out_path) << cons_newick << "\n";
        if (resolve) {
            long polytomies = 0;     // unrooted: a root with three children is not one
            for (int v = 0; v < forest.nodes(); ++v)
                polytomies += forest.child_off[v + 1] - forest.child_off[v] > (forest.parent[v] < 0 ? 3 : 2);
            const int before = forest.nodes();
            forest = resolve_polytomies(forest);
            log("resolved " + std::to_string(polytomies) + " polytomies arbitrarily (" + std::to_string(forest.nodes() - before) +
                " nodes added, counting the split of three-child roots)");
        }
        ConsensusPartitioner cons(parse_newick(cons_newick), taxa);
        cons.set_objective(opt.kappa, opt.power);
        cons.set_init_kappa(opt.init_kappa);
        std::unique_ptr<Trace> trace;
        if (!trace_path.empty()) {
            trace = std::make_unique<Trace>(trace_path);
            for (int t = 0; t < taxa.size(); ++t) trace->line("T " + std::to_string(t) + " " + taxa.label[t]);
        }

        std::unique_ptr<std::FILE, int (*)(std::FILE*)> sublog(nullptr, std::fclose);
        if (!sublog_path.empty()) {
            sublog.reset(std::fopen(sublog_path.c_str(), "w"));
            if (!sublog) throw std::runtime_error("cannot write " + sublog_path);
            opt.subproblem_log = sublog.get();
        }

        log("searching with " + std::to_string(opt.threads) + " thread(s)");
        Search search(forest, cons, trace.get(), replay.get(), opt);
        const SpeciesTree tree = search.run();
        sublog.reset();
        std::ofstream(out_path) << tree.newick(taxa) << "\n";
        log("search finished; species tree written to " + out_path);

        const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        const SearchStats st = search.stats();
        std::fprintf(stderr, "taxa %d, gene trees %d, nodes %d; subproblems %ld, FM passes %ld, moves %ld; %d threads; %.2f s\n",
                     taxa.size(), forest.k, forest.nodes(), st.subproblems, st.passes, st.moves, st.threads, secs);
        std::fprintf(stderr, "time (summed over threads): initial bipartitions %.2f s, FM passes %.2f s, subproblem setup %.2f s\n",
                     st.t_initial, st.t_fm, st.t_setup);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
