#include "newick.hpp"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string_view>

namespace wqfm {
namespace {

constexpr std::string_view kDelims = ",():;[";

class Parser {
public:
    explicit Parser(const std::string& s) : s_(s) {}

    ParsedTree run() {
        std::vector<std::vector<int>> open;   // children of each currently open '('
        for (;;) {
            skip();
            if (i_ >= s_.size()) fail("missing ')'");
            const char c = s_[i_];
            if (c == '(') {
                open.emplace_back();
                ++i_;
            } else if (c == ',') {
                if (open.empty()) fail("',' outside parentheses");
                ++i_;
            } else if (c == ')') {
                if (open.empty()) fail("unbalanced ')'");
                ++i_;
                std::vector<int> kids = std::move(open.back());
                open.pop_back();
                if (kids.empty()) fail("empty subtree");
                const int v = add_node("");
                t_.children[v].assign(kids.rbegin(), kids.rend());
                for (int k : kids) t_.parent[k] = v;
                t_.support[v] = numeric(label());  // an internal label that is a number: the branch's support
                t_.length[v] = branch_length();
                if (open.empty()) {
                    t_.root = v;
                    skip();
                    if (i_ < s_.size() && s_[i_] == ';') ++i_;
                    skip();
                    if (i_ != s_.size()) fail("text after the end of the tree");
                    return std::move(t_);
                }
                open.back().push_back(v);
                expect_separator();
            } else if (c == ';') {
                fail("';' before the tree is closed");
            } else {
                std::string lab = label();
                if (lab.empty()) fail("empty label");
                if (open.empty()) fail("leaf outside parentheses");
                const int v = add_node(std::move(lab));
                t_.length[v] = branch_length();
                open.back().push_back(v);
                expect_separator();
            }
        }
    }

private:
    [[noreturn]] void fail(const char* what) const {
        throw std::runtime_error(std::string("newick: ") + what + " at position " + std::to_string(i_));
    }

    void skip() {
        for (;;) {
            while (i_ < s_.size() && std::isspace(static_cast<unsigned char>(s_[i_]))) ++i_;
            if (i_ < s_.size() && s_[i_] == '[') {
                const size_t j = s_.find(']', i_);
                if (j == std::string::npos) fail("unterminated comment");
                i_ = j + 1;
            } else {
                return;
            }
        }
    }

    bool plain(char c) const {
        return kDelims.find(c) == std::string_view::npos && !std::isspace(static_cast<unsigned char>(c));
    }

    std::string label() {
        skip();
        std::string out;
        if (i_ < s_.size() && s_[i_] == '\'') {
            ++i_;
            for (;;) {
                if (i_ >= s_.size()) fail("unterminated quoted label");
                if (s_[i_] == '\'') {
                    if (i_ + 1 < s_.size() && s_[i_ + 1] == '\'') { out += '\''; i_ += 2; continue; }
                    ++i_;
                    return out;
                }
                out += s_[i_++];
            }
        }
        while (i_ < s_.size() && plain(s_[i_])) out += s_[i_++];
        return out;
    }

    double branch_length() {                 // NaN if absent or not a number
        skip();
        double out = std::nan("");
        if (i_ < s_.size() && s_[i_] == ':') {
            ++i_;
            skip();
            std::string tok;
            while (i_ < s_.size() && plain(s_[i_])) tok += s_[i_++];
            out = numeric(tok);
        }
        skip();
        return out;
    }

    static double numeric(const std::string& tok) {   // the whole token as a number, else NaN
        if (tok.empty()) return std::nan("");
        char* end = nullptr;
        const double x = std::strtod(tok.c_str(), &end);
        return end == tok.c_str() + tok.size() ? x : std::nan("");
    }

    // after a complete subtree the next token must be ',' or ')'; catches unquoted labels with spaces
    void expect_separator() {
        skip();
        if (i_ >= s_.size() || (s_[i_] != ',' && s_[i_] != ')')) fail("expected ',' or ')'");
    }

    int add_node(std::string lab) {
        const int v = static_cast<int>(t_.parent.size());
        t_.parent.push_back(-1);
        t_.children.emplace_back();
        t_.label.push_back(std::move(lab));
        t_.support.push_back(std::nan(""));
        t_.length.push_back(std::nan(""));
        return v;
    }

    const std::string& s_;
    size_t i_ = 0;
    ParsedTree t_;
};

}  // namespace

ParsedTree parse_newick(const std::string& text) { return Parser(text).run(); }

}  // namespace wqfm
