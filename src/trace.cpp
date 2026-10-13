#include "trace.hpp"

#include <fstream>
#include <stdexcept>

namespace wqfm {

Trace::Trace(const std::string& path) : f_(std::fopen(path.c_str(), "w")) {
    if (!f_) throw std::runtime_error("cannot write trace file " + path);
}

Trace::~Trace() { std::fclose(f_); }

void Trace::line(const std::string& s) {
    std::fputs(s.c_str(), f_);
    std::fputc('\n', f_);
}

std::string Trace::num(double x) {
    char buf[40];
    std::snprintf(buf, sizeof buf, "%.17g", x);
    return buf;
}

std::pair<std::vector<int8_t>, std::vector<int8_t>> parse_bits(const std::string& bits) {
    std::pair<std::vector<int8_t>, std::vector<int8_t>> out;
    bool dummies = false;
    for (char c : bits) {
        if (c == '|') { dummies = true; continue; }
        if (c != '0' && c != '1') throw std::runtime_error("bad partition bits: " + bits);
        (dummies ? out.second : out.first).push_back(static_cast<int8_t>(c - '0'));
    }
    return out;
}

Replay Replay::load(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read replay trace " + path);
    Replay r;
    Part* part = nullptr;
    Sub* sub = nullptr;
    bool pass_open = false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.size() < 2 || line[1] != ' ') continue;
        const std::string rest = line.substr(2);
        switch (line[0]) {
            case 'T': {
                const size_t sp = rest.find(' ');
                r.taxa.emplace_back(std::stoi(rest.substr(0, sp)), rest.substr(sp + 1));
                break;
            }
            case 'P': part = &r.parts[rest]; *part = Part{}; break;
            case 'C': if (part) part->cand.push_back(std::stod(rest)); break;
            case 'R': if (part) part->result = rest; break;
            case 'S': sub = &r.subs[rest]; *sub = Sub{}; pass_open = false; break;
            case 'I': if (sub) sub->init = rest; break;
            case 'M': {
                if (!sub) break;
                if (!pass_open) { sub->passes.emplace_back(); pass_open = true; }
                Move m;
                size_t pos = 0;
                m.dummy = std::stoi(rest, &pos);
                std::string tail = rest.substr(pos);
                m.index = std::stoi(tail, &pos); tail = tail.substr(pos);
                m.gain = std::stod(tail, &pos); tail = tail.substr(pos);
                m.cg = std::stod(tail);
                sub->passes.back().moves.push_back(m);
                break;
            }
            case 'E':
                if (!sub) break;
                if (!pass_open) sub->passes.emplace_back();
                sub->passes.back().max_cg = std::stoi(rest);
                pass_open = false;
                break;
            case 'F': if (sub) sub->fin = rest; break;
            default: break;
        }
    }
    return r;
}

}  // namespace wqfm
