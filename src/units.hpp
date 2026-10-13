#pragma once

#include <cstdint>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace wqfm {

// A dummy taxon stands for one side of an earlier bipartition: its children are that side's real
// taxa and (nested) dummy taxa. `flat` lists every real taxon under it.
struct Dummy {
    std::vector<int> reals;
    std::vector<int> subs;     // DummyPool ids
    std::vector<int> flat;
};

// Thread-safe for concurrent add() and reads once reserve() has fixed the capacity (entries never move).
class DummyPool {
public:
    void reserve(size_t n) {
        std::lock_guard<std::mutex> lk(m_);
        pool_.reserve(n);
        fixed_ = true;
    }

    int add(std::vector<int> reals, std::vector<int> subs) {
        std::lock_guard<std::mutex> lk(m_);
        if (fixed_ && pool_.size() == pool_.capacity()) throw std::logic_error("DummyPool: reserved capacity exceeded");
        Dummy d;
        d.flat = reals;
        for (int s : subs) d.flat.insert(d.flat.end(), pool_[s].flat.begin(), pool_[s].flat.end());
        d.reals = std::move(reals);
        d.subs = std::move(subs);
        pool_.push_back(std::move(d));
        return static_cast<int>(pool_.size()) - 1;
    }

    const Dummy& operator[](int id) const { return pool_[id]; }

    // Java DummyTaxon.calcDivCoeffs (nested normalization, Eq. 1 of the paper): the coefficient of a
    // real taxon is the product of the child counts on the path to it; its weight is 1/coefficient.
    void coeffs(int id, double mult, std::vector<double>& coeff) const {
        const Dummy& d = pool_[id];
        const double sz = static_cast<double>(d.reals.size() + d.subs.size());
        for (int r : d.reals) coeff[r] = sz * mult;
        for (int s : d.subs) coeffs(s, sz * mult, coeff);
    }

    // "{r1,r2|{...};{...}}" - shared with the Java trace
    std::string desc(int id) const {
        const Dummy& d = pool_[id];
        std::string s = "{";
        for (size_t i = 0; i < d.reals.size(); ++i) s += (i ? "," : "") + std::to_string(d.reals[i]);
        s += "|";
        for (size_t i = 0; i < d.subs.size(); ++i) s += (i ? ";" : "") + desc(d.subs[i]);
        return s + "}";
    }

private:
    std::vector<Dummy> pool_;
    std::mutex m_;
    bool fixed_ = false;
};

// The units of one divide step, in the Java order (it determines tie-breaking), and their bipartition.
struct Subproblem {
    std::vector<int> reals;          // taxon ids
    std::vector<int> dummies;        // DummyPool ids
    std::vector<int8_t> part_r;      // partition (0/1) of each real taxon
    std::vector<int8_t> part_d;      // partition of each dummy taxon

    int size() const { return static_cast<int>(reals.size() + dummies.size()); }
};

// trace keys, shared with the Java trace
inline std::string sub_key(const Subproblem& sp, const DummyPool& pool) {
    std::string s = "R:";
    for (size_t i = 0; i < sp.reals.size(); ++i) s += (i ? "," : "") + std::to_string(sp.reals[i]);
    s += " D:";
    for (size_t i = 0; i < sp.dummies.size(); ++i) s += (i ? ";" : "") + pool.desc(sp.dummies[i]);
    return s;
}

inline std::string part_bits(const std::vector<int8_t>& r, const std::vector<int8_t>& d) {
    std::string s;
    for (int8_t x : r) s += static_cast<char>('0' + x);
    s += '|';
    for (int8_t x : d) s += static_cast<char>('0' + x);
    return s;
}

}  // namespace wqfm
