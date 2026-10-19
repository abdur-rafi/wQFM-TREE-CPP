#include "wqfm.hpp"

#include "wscorer.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace wqfm {

int SpeciesTree::add(int t, int d) {
    adj.emplace_back();
    taxon.push_back(t);
    dummy.push_back(d);
    return static_cast<int>(adj.size()) - 1;
}

std::string SpeciesTree::newick(const Taxa& taxa) const {
    auto is_leaf = [&](int v) { return taxon[v] >= 0 || dummy[v] >= 0; };
    int root = 0;
    while (root < static_cast<int>(adj.size()) && is_leaf(root)) ++root;
    if (root == static_cast<int>(adj.size())) throw std::logic_error("species tree has no internal node");

    struct Frame { int v, parent; size_t next; };
    std::string out = "(";
    std::vector<Frame> st{{root, -1, 0}};
    while (!st.empty()) {
        Frame& f = st.back();
        while (f.next < adj[f.v].size() && adj[f.v][f.next] == f.parent) ++f.next;
        if (f.next == adj[f.v].size()) {
            out += ')';
            st.pop_back();
            continue;
        }
        const int c = adj[f.v][f.next++], from = f.v;
        if (out.back() != '(') out += ',';
        if (is_leaf(c)) {
            out += taxon[c] >= 0 ? taxa.label[taxon[c]] : "#dummy" + std::to_string(dummy[c]);
        } else {
            out += '(';
            st.push_back({c, from, 0});   // invalidates f
        }
    }
    return out + ";";
}

namespace {

using Clock = std::chrono::steady_clock;

double since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

SpeciesTree star(const Subproblem& sp) {
    SpeciesTree t;
    const int c = t.add(-1, -1);
    auto leaf = [&](int taxon, int dummy) {
        const int v = t.add(taxon, dummy);
        t.adj[c].push_back(v);
        t.adj[v].push_back(c);
    };
    for (int r : sp.reals) leaf(r, -1);
    for (int d : sp.dummies) leaf(-1, d);
    return t;
}

// Replace dummy leaf x0 of t0 and dummy leaf x1 of t1 by one edge between their neighbours.
SpeciesTree join(const SpeciesTree& t0, int x0, const SpeciesTree& t1, int x1) {
    auto find = [](const SpeciesTree& t, int x) {
        for (size_t v = 0; v < t.adj.size(); ++v)
            if (t.dummy[v] == x) return static_cast<int>(v);
        throw std::logic_error("dummy leaf not found while joining subtrees");
    };
    const int p0 = find(t0, x0), p1 = find(t1, x1);
    if (t0.adj[p0].size() != 1 || t1.adj[p1].size() != 1) throw std::logic_error("dummy taxon is not a leaf");
    SpeciesTree out;
    std::vector<int> m0(t0.adj.size(), -1), m1(t1.adj.size(), -1);
    for (size_t v = 0; v < t1.adj.size(); ++v) if (static_cast<int>(v) != p1) m1[v] = out.add(t1.taxon[v], t1.dummy[v]);
    for (size_t v = 0; v < t0.adj.size(); ++v) if (static_cast<int>(v) != p0) m0[v] = out.add(t0.taxon[v], t0.dummy[v]);
    for (size_t v = 0; v < t1.adj.size(); ++v)
        if (m1[v] >= 0) for (int w : t1.adj[v]) if (w != p1) out.adj[m1[v]].push_back(m1[w]);
    for (size_t v = 0; v < t0.adj.size(); ++v)
        if (m0[v] >= 0) for (int w : t0.adj[v]) if (w != p0) out.adj[m0[v]].push_back(m0[w]);
    const int q0 = m0[t0.adj[p0][0]], q1 = m1[t1.adj[p1][0]];
    out.adj[q1].push_back(q0);
    out.adj[q0].push_back(q1);
    return out;
}

}  // namespace

Search::Search(const Forest& f, const ConsensusPartitioner& cons, Trace* trace, const Replay* replay,
               const SearchOptions& opt)
    : F_(f), cons_(cons), trace_(trace), replay_(replay), compress_(opt.compress), relative_eps_(opt.relative_eps),
      progress_secs_(opt.progress_secs), t_start_(std::chrono::steady_clock::now()), sublog_(opt.subproblem_log),
      kappa_(opt.kappa), power_(opt.power), gamma_(opt.dummy_gamma), weighting_(opt.weighting), smin_(opt.support_min), smax_(opt.support_max) {
    if (sublog_) std::fprintf(sublog_, "# reals\tdummies\tinitial_score\tfinal_score\tfm_passes\tchanged\n");
    next_report_ns_ = static_cast<long long>(opt.progress_secs * 1e9);
    if (opt.threads > 1) threads_ = std::make_unique<ThreadPool>(opt.threads);
    tasks_ = opt.threads > 1 && !trace && !replay;
    pool_.reserve(2 * static_cast<size_t>(f.n) + 16);   // two dummies per divide step, at most n steps
}

SearchStats Search::stats() const {
    std::lock_guard<std::mutex> lk(mu_);
    SearchStats st;
    st.subproblems = subproblems_;
    st.passes = passes_;
    st.moves = moves_;
    st.t_setup = t_setup_;
    st.t_initial = t_initial_;
    st.t_fm = t_fm_;
    st.threads = threads_ ? threads_->threads() : 1;
    return st;
}

void Search::add_time(double& total, double secs) {
    std::lock_guard<std::mutex> lk(mu_);
    total += secs;
}

// A fresh scorer per subproblem, freed when the subproblem is done. Reusing them saved no measurable
// time, but a reused scorer kept the buffers of the largest forest it had ever held, so with several
// threads peak memory grew by up to 60% and depended on scheduling.
std::unique_ptr<IScorer> Search::acquire() {
    std::unique_ptr<IScorer> sc;
    if (weighting_) {
        WScorer::Options o;
        o.support = weighting_ & 1;
        o.length = weighting_ & 2;
        o.smin = smin_;
        o.smax = smax_;
        sc = std::make_unique<WScorer>(F_, compress_, o);
    } else {
        sc = std::make_unique<Scorer>(F_, compress_);
    }
    sc->set_kappa(power_ == 0 ? kappa_ : 1.0);   // a normalised objective combines kappa 0 and 1 itself
    sc->set_dummy_gamma(gamma_);
    return sc;
}

double Search::objective_value(IScorer& sc) {
    if (power_ == 0) return sc.score();
    sc.set_kappa(0);
    const double s0 = sc.score();
    sc.set_kappa(1);
    const double T = s0 - sc.score();
    return T > 0 ? (s0 - kappa_ * T) / std::pow(T, power_) : 0.0;
}

// The objective (2F - kappa T) / T^power and every unit's gain. With power 0 the scorer computes it
// directly; otherwise two passes give 2F (kappa 0) and 2F - T (kappa 1), hence F, T and their changes.
double Search::objective_gains(IScorer& sc, std::vector<double>& gr, std::vector<double>& gd) {
    if (power_ == 0) return sc.gains(gr, gd);
    std::vector<double> fr, fd;
    sc.set_kappa(0);
    const double s0 = sc.gains(fr, fd);
    sc.set_kappa(1);
    const double s1 = sc.gains(gr, gd);
    auto obj = [&](double f2, double t) { return t > 0 ? (f2 - kappa_ * t) / std::pow(t, power_) : 0.0; };
    const double T = s0 - s1, cur = obj(s0, T);
    for (size_t i = 0; i < gr.size(); ++i) gr[i] = obj(s0 + fr[i], T + fr[i] - gr[i]) - cur;
    for (size_t j = 0; j < gd.size(); ++j) gd[j] = obj(s0 + fd[j], T + fd[j] - gd[j]) - cur;
    return cur;
}

void Search::release(std::unique_ptr<IScorer> sc) {
    sc.reset();
}

SpeciesTree Search::run() {
    Subproblem top;
    for (int t = 0; t < F_.n; ++t) top.reals.push_back(t);
    std::unique_ptr<IScorer> sc = prepare(top, nullptr, 0);
    return recurse(top, std::move(sc));
}

std::unique_ptr<IScorer> Search::prepare(Subproblem& sp, const IScorer* parent, int side) {
    std::unique_ptr<IScorer> sc = acquire();
    sc->set_pool(threads_.get());
    auto t0 = Clock::now();
    if (parent && compress_) sc->setup_child(sp, *parent, side);
    else sc->setup(sp, pool_);
    add_time(t_setup_, since(t0));
    t0 = Clock::now();
    cons_.make(sp, true, *sc, pool_, trace_, replay_);   // Java passes allowSingleton = true
    add_time(t_initial_, since(t0));
    return sc;
}

// A subproblem with more than 3 units is split in two, so a binary tree on n taxa has n - 3 of them:
// the count is a good measure of progress (the large ones come first, so time is not linear in it).
void Search::report_progress() {
    if (progress_secs_ <= 0) return;
    const long long now = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t_start_).count();
    long long due = next_report_ns_.load();
    if (now < due || !next_report_ns_.compare_exchange_strong(due, now + static_cast<long long>(progress_secs_ * 1e9))) return;
    const long total = std::max(1, F_.n - 3);
    const long done = std::min(subproblems_.load(), total);
    std::fprintf(stderr, "[%7.1f s] search: %ld of ~%ld subproblems (%.0f%%), %ld FM passes, %ld moves\n",
                 now / 1e9, done, total, 100.0 * done / total, passes_.load(), moves_.load());
}

SpeciesTree Search::recurse(Subproblem& sp, std::unique_ptr<IScorer> sc) {
    ++subproblems_;
    report_progress();
    const Replay::Sub* rs = nullptr;
    if (trace_ || replay_) {
        const std::string key = sub_key(sp, pool_);
        if (trace_) {
            trace_->line("S " + key);
            trace_->line("I " + part_bits(sp.part_r, sp.part_d));
        }
        if (replay_) {
            auto it = replay_->subs.find(key);
            if (it != replay_->subs.end()) rs = &it->second;
        }
    }

    const auto t0 = Clock::now();
    sc->set_partition(sp.part_r, sp.part_d);
    double initial_score = 0;
    std::vector<int8_t> init_r, init_d;
    if (sublog_) {
        initial_score = objective_value(*sc);
        init_r = sp.part_r;
        init_d = sp.part_d;
    }
    size_t pass = 0;
    for (int itr = 0;;) {                      // QFM2.recurse: while(oneInteration) if(++itr > 5) break
        const Replay::Pass* rp = nullptr;
        if (rs) {
            if (pass >= rs->passes.size()) throw std::runtime_error("replay: more FM passes than in the Java trace");
            rp = &rs->passes[pass];
        }
        ++pass;
        if (!fm_pass(sp, *sc, rp)) break;
        if (++itr > 5) break;
    }
    add_time(t_fm_, since(t0));
    for (size_t i = 0; i < sp.reals.size(); ++i) sp.part_r[i] = static_cast<int8_t>(sc->part_real(static_cast<int>(i)));
    for (size_t j = 0; j < sp.dummies.size(); ++j) sp.part_d[j] = static_cast<int8_t>(sc->part_dummy(static_cast<int>(j)));
    if (trace_) trace_->line("F " + part_bits(sp.part_r, sp.part_d));
    if (sublog_) {
        const double final_score = objective_value(*sc);
        bool same = true, flipped = true;         // the same bipartition, or its sides swapped
        for (size_t i = 0; i < init_r.size(); ++i) { same &= init_r[i] == sp.part_r[i]; flipped &= init_r[i] != sp.part_r[i]; }
        for (size_t j = 0; j < init_d.size(); ++j) { same &= init_d[j] == sp.part_d[j]; flipped &= init_d[j] != sp.part_d[j]; }
        std::lock_guard<std::mutex> lk(log_mu_);
        std::fprintf(sublog_, "%zu\t%zu\t%.17g\t%.17g\t%zu\t%d\n", sp.reals.size(), sp.dummies.size(), initial_score,
                     final_score, pass, same || flipped ? 0 : 1);
    }

    // divide: each side gets a new dummy taxon standing for the other side
    std::vector<int> rts[2], dts[2];
    for (size_t i = 0; i < sp.reals.size(); ++i) rts[sp.part_r[i]].push_back(sp.reals[i]);
    for (size_t j = 0; j < sp.dummies.size(); ++j) dts[sp.part_d[j]].push_back(sp.dummies[j]);
    Subproblem child[2];
    std::unique_ptr<IScorer> child_sc[2];
    int fresh[2];
    for (int s = 0; s < 2; ++s) {
        fresh[s] = pool_.add(rts[1 - s], dts[1 - s]);
        child[s].reals = rts[s];
        child[s].dummies = dts[s];
        child[s].dummies.push_back(fresh[s]);
        if (child[s].size() > 3) child_sc[s] = prepare(child[s], sc.get(), s);
    }
    release(std::move(sc));                      // done with this subproblem's scorer

    // fork-join: child 1 as a task (picked up by an idle thread), child 0 here
    SpeciesTree t[2];
    std::shared_ptr<ThreadPool::Task> task;
    if (tasks_ && child[0].size() > 3 && child[1].size() > 3)
        task = threads_->submit([this, &child, &child_sc, &t] { t[1] = recurse(child[1], std::move(child_sc[1])); });
    try {
        t[0] = child[0].size() > 3 ? recurse(child[0], std::move(child_sc[0])) : star(child[0]);
    } catch (...) {
        if (task) {                              // the task refers to this frame: let it finish first
            try { threads_->wait(task); } catch (...) {}
        }
        throw;
    }
    if (task) threads_->wait(task);
    else t[1] = child[1].size() > 3 ? recurse(child[1], std::move(child_sc[1])) : star(child[1]);
    return join(t[0], fresh[0], t[1], fresh[1]);
}

// One FM pass (QFM2.oneInteration + swapMax): repeatedly move the unlocked unit with the largest gain
// out of a side with more than 2 units, then keep the prefix with the largest cumulative gain.
bool Search::fm_pass(Subproblem& sp, IScorer& sc, const Replay::Pass* rp) {
    const int R = static_cast<int>(sp.reals.size()), D = static_cast<int>(sp.dummies.size());
    int cnt[2] = {0, 0};
    for (int i = 0; i < R; ++i) ++cnt[sc.part_real(i)];
    for (int j = 0; j < D; ++j) ++cnt[sc.part_dummy(j)];
    const bool singleton = cnt[0] == 1 || cnt[1] == 1;
    std::vector<char> lock_r(R, 0), lock_d(D, 0);
    struct Move { bool dummy; int index; };
    std::vector<Move> swaps;
    std::vector<double> gr, gd;
    double cg = 0, max_cg = 0, moved = 0;   // moved: sum of |gain| of the moves so far
    int max_idx = -1;

    for (size_t step = 0;; ++step) {
        objective_gains(sc, gr, gd);
        int best = -1;
        double bg = 0;
        bool dum = false;
        for (int i = 0; i < R; ++i) {
            if (lock_r[i] || cnt[sc.part_real(i)] <= 2) continue;
            if (best == -1 || bg < gr[i]) { best = i; bg = gr[i]; }
        }
        for (int j = 0; j < D; ++j) {
            if (lock_d[j] || cnt[sc.part_dummy(j)] <= 2) continue;
            if (best == -1 || bg < gd[j]) { best = j; bg = gd[j]; dum = true; }
        }
        const int own = best;
        const bool own_dum = dum;
        const double own_gain = bg;
        if (rp) {                              // follow Java's move, keep ours for the trace
            if (step < rp->moves.size()) {
                const Replay::Move& m = rp->moves[step];
                dum = m.dummy != 0;
                best = m.index;
                if (best < 0 || best >= (dum ? D : R) || (dum ? lock_d[best] : lock_r[best]))
                    throw std::runtime_error("replay: invalid move in the Java trace");
                bg = dum ? gd[best] : gr[best];
            } else {
                best = -1;
            }
        }
        if (best == -1) break;

        const int p = dum ? sc.part_dummy(best) : sc.part_real(best);
        if (dum) { sc.swap_dummy(best); lock_d[best] = 1; }
        else { sc.swap_real(best); lock_r[best] = 1; }
        --cnt[p];
        ++cnt[1 - p];
        swaps.push_back({dum, best});
        cg += bg;
        moved += std::fabs(bg);
        ++moves_;
        if (trace_) {
            std::string s = "M " + std::to_string(dum ? 1 : 0) + " " + std::to_string(best) + " " + Trace::num(bg) + " " + Trace::num(cg);
            if (rp) s += " " + std::to_string(own_dum ? 1 : 0) + " " + std::to_string(own) + " " + Trace::num(own_gain);
            trace_->line(s);
        }
        if (singleton && max_idx == -1) { max_cg = cg; max_idx = static_cast<int>(swaps.size()) - 1; }
        const double eps = relative_eps_ ? 1e-5 + 1e-12 * moved : 1e-5;
        if (cg > max_cg && std::fabs(max_cg - cg) > eps) { max_cg = cg; max_idx = static_cast<int>(swaps.size()) - 1; }
    }
    if (rp) max_idx = rp->max_cg;
    if (trace_) trace_->line("E " + std::to_string(max_idx));
    ++passes_;

    auto undo = [&](const Move& m) { m.dummy ? sc.swap_dummy(m.index) : sc.swap_real(m.index); };
    if (max_idx == -1) {
        if (static_cast<int>(swaps.size()) != R + D)
            for (int i = static_cast<int>(swaps.size()) - 1; i >= 0; --i) undo(swaps[i]);
        return false;
    }
    for (int i = static_cast<int>(swaps.size()) - 1; i > max_idx; --i) undo(swaps[i]);
    return true;
}

}  // namespace wqfm
