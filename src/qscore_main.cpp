// qscore: quartet score of species trees against gene trees (ASTRAL's quartet score, counted exactly).
//
// usage: qscore -g <gene trees> [-t threads] <species tree file>...
//
// For each species tree file (its first tree) prints one tab-separated line: the gene-tree quartets
// the species tree agrees with, the quartets the gene trees resolve, the normalised score (their
// ratio) and the file name. Trees are unrooted; gene trees may have polytomies and missing taxa.

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "forest.hpp"
#include "newick.hpp"
#include "qscore.hpp"

using namespace wqfm;

namespace {

std::vector<ParsedTree> read_trees(const std::string& path, bool first_only) {
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
        if (first_only) break;
    }
    if (trees.empty()) throw std::runtime_error("no trees in " + path);
    return trees;
}

void add_taxa(const ParsedTree& t, Taxa& taxa) {
    for (size_t v = 0; v < t.label.size(); ++v) if (t.children[v].empty()) taxa.add(t.label[v]);
}

[[noreturn]] void usage() {
    std::fprintf(stderr, "usage: qscore -g <gene trees> [-t threads] <species tree file>...\n");
    std::exit(2);
}

}  // namespace

int main(int argc, char** argv) {
    std::string genes_path;
    std::vector<std::string> species_paths;
    int threads = 1;
    for (int a = 1; a < argc; ++a) {
        const std::string arg = argv[a];
        auto value = [&]() -> std::string { if (a + 1 >= argc) usage(); return argv[++a]; };
        if (arg == "-g") genes_path = value();
        else if (arg == "-t") threads = std::max(1, std::atoi(value().c_str()));
        else if (!arg.empty() && arg[0] == '-') usage();
        else species_paths.push_back(arg);
    }
    if (genes_path.empty() || species_paths.empty()) usage();
    try {
        const std::vector<ParsedTree> gene_trees = read_trees(genes_path, false);
        std::vector<ParsedTree> species;
        for (const auto& p : species_paths) species.push_back(read_trees(p, true)[0]);
        Taxa taxa;
        for (const auto& t : gene_trees) add_taxa(t, taxa);
        for (const auto& t : species) add_taxa(t, taxa);
        const Forest genes = make_forest(gene_trees, taxa);
        for (size_t s = 0; s < species.size(); ++s) {
            const Forest sp = make_forest({species[s]}, taxa);
            QuartetCounts q;
            try {
                q = quartet_support(genes, sp, threads);
            } catch (const std::exception& e) {
                throw std::runtime_error(species_paths[s] + ": " + e.what());
            }
            std::printf("%lld\t%lld\t%.10f\t%s\n", static_cast<long long>(q.agree), static_cast<long long>(q.resolved),
                        q.resolved ? static_cast<double>(q.agree) / static_cast<double>(q.resolved) : 0.0,
                        species_paths[s].c_str());
            std::fflush(stdout);
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
