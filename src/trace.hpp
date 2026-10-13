#pragma once

#include <cstdint>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace wqfm {

// Step-by-step trace of the search, in the format written by the instrumented Java build
// (tools/java/Trace.java), so the two can be compared line by line:
//   T id label               taxon order
//   P key / C score / R bits initial-bipartition candidates of a subproblem, and the chosen one
//   S key / I bits           subproblem start and its initial bipartition
//   M d i gain cg [od oi og] FM move: dummy?, unit index, gain, cumulative gain [replay: own choice]
//   E index                  end of an FM pass: prefix kept (-1 = none)
//   F bits                   final bipartition
// key = "R:<real ids> D:<dummy descriptors>", bits = "<real sides>|<dummy sides>".
class Trace {
public:
    explicit Trace(const std::string& path);
    ~Trace();
    Trace(const Trace&) = delete;
    Trace& operator=(const Trace&) = delete;

    void line(const std::string& s);
    static std::string num(double x);    // 17 significant digits: round-trips exactly

private:
    FILE* f_;
};

// A Java trace, for replay: the C++ search then follows Java's decisions (initial bipartitions, FM
// moves, kept prefixes) while computing and tracing its own values, so every number can be compared
// even where floating-point rounding would break a near-tie differently.
struct Replay {
    struct Move { int dummy = 0, index = 0; double gain = 0, cg = 0; };
    struct Pass { std::vector<Move> moves; int max_cg = -1; };
    struct Sub { std::string init, fin; std::vector<Pass> passes; };
    struct Part { std::vector<double> cand; std::string result; };

    std::unordered_map<std::string, Sub> subs;
    std::unordered_map<std::string, Part> parts;
    std::vector<std::pair<int, std::string>> taxa;   // T lines

    static Replay load(const std::string& path);
};

std::pair<std::vector<int8_t>, std::vector<int8_t>> parse_bits(const std::string& bits);

}  // namespace wqfm
