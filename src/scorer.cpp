#include "scorer.hpp"

#include <algorithm>
#include <stdexcept>

namespace wqfm {

// ---------------------------------------------------------------- PolyEval
//
// With s_m = +1 for m = 0 and -1 for dummy components:
//   EA = sum_{i<j} PA_ij = sum_m s_m (TA_m^2 - QA_m)/2,  EB likewise
//   S  = sum_{i<j} PA_ij (P - PB_i - PB_j)               = P EA - sum_m s_m (TA_m R_m - R2_m)
//   U  = sum_{i<j} PA_ij sum_{k<l, k,l not in {i,j}} PB_kl = EA EB - sum_i SB_i g_i + sum_{i<j} PA_ij PB_ij
// where g_i = sum_{j!=i} PA_ij and SB_i = sum_{l!=i} PB_il, which expand into the M1..M4 power sums.

namespace {
inline double sgn(int m) { return m == 0 ? 1.0 : -1.0; }
}

void PolyEval::resize(int d, int ma, int mb) {
    D = d; mA = ma; mB = mb;
    al.assign(static_cast<size_t>(ma) * d, 0.0);
    be.assign(static_cast<size_t>(mb) * d, 0.0);
    rho.assign(d, 0.0);
    agg.assign(4 * ma + 2 * mb + 4 * ma * mb + 1, 0.0);
}

void PolyEval::build() {
    std::fill(agg.begin(), agg.end(), 0.0);
    for (int i = 0; i < D; ++i) add(i, 1.0);
}

void PolyEval::add(int i, double s) {
    double* TA = agg.data();
    double* QA = TA + mA;
    double* R = QA + mA;
    double* R2 = R + mA;
    double* TB = R2 + mA;
    double* QB = TB + mB;
    double* M1 = QB + mB;
    double* M2 = M1 + mA * mB;
    double* M3 = M2 + mA * mB;
    double* M4 = M3 + mA * mB;
    double* P = M4 + mA * mB;

    double pb = -rho[i];
    for (int n = 0; n < mB; ++n) {
        const double b = be[static_cast<size_t>(n) * D + i];
        pb += sgn(n) * b * b;
    }
    pb *= 0.5;
    for (int m = 0; m < mA; ++m) {
        const double a = al[static_cast<size_t>(m) * D + i];
        if (a == 0) continue;
        TA[m] += s * a;
        QA[m] += s * a * a;
        R[m] += s * pb * a;
        R2[m] += s * pb * a * a;
        for (int n = 0; n < mB; ++n) {
            const double b = be[static_cast<size_t>(n) * D + i];
            if (b == 0) continue;
            const int k = m * mB + n;
            M1[k] += s * a * b;
            M2[k] += s * a * a * b;
            M3[k] += s * a * b * b;
            M4[k] += s * a * a * b * b;
        }
    }
    for (int n = 0; n < mB; ++n) {
        const double b = be[static_cast<size_t>(n) * D + i];
        TB[n] += s * b;
        QB[n] += s * b * b;
    }
    *P += s * pb;
}

double PolyEval::eval() const {
    const double* TA = agg.data();
    const double* QA = TA + mA;
    const double* R = QA + mA;
    const double* R2 = R + mA;
    const double* TB = R2 + mA;
    const double* QB = TB + mB;
    const double* M1 = QB + mB;
    const double* M2 = M1 + mA * mB;
    const double* M3 = M2 + mA * mB;
    const double* M4 = M3 + mA * mB;
    const double P = *(M4 + mA * mB);

    double EA = 0, EB = 0, sr = 0, X = 0;
    for (int m = 0; m < mA; ++m) {
        EA += sgn(m) * (TA[m] * TA[m] - QA[m]) / 2;
        sr += sgn(m) * (TA[m] * R[m] - R2[m]);
    }
    for (int n = 0; n < mB; ++n) EB += sgn(n) * (TB[n] * TB[n] - QB[n]) / 2;
    for (int m = 0; m < mA; ++m)
        for (int n = 0; n < mB; ++n) {
            const int k = m * mB + n;
            X += sgn(m) * sgn(n) *
                 (TB[n] * TA[m] * M1[k] - TB[n] * M2[k] - TA[m] * M3[k] + M4[k] - (M1[k] * M1[k] - M4[k]) / 2);
        }
    const double S = P * EA - sr;
    const double U = EA * EB - X;
    return S + U / 2;
}

// ---------------------------------------------------------------- Scorer

void Scorer::setup(const Subproblem& sp, const DummyPool& pool) {
    const int n = F_.n, N = F_.nodes();
    d_ = static_cast<int>(sp.dummies.size());
    reals_ = sp.reals;
    unit_r_.assign(n, -1);
    unit_d_.assign(n, -1);
    std::vector<double> coeff(n, 0.0);
    for (size_t i = 0; i < sp.reals.size(); ++i) {
        unit_r_[sp.reals[i]] = static_cast<int>(i);
        coeff[sp.reals[i]] = 1.0;
    }
    for (int j = 0; j < d_; ++j) {
        pool.coeffs(sp.dummies[j], 1.0, coeff);
        for (int t : pool[sp.dummies[j]].flat) unit_d_[t] = j;
    }
    w_.assign(n, 0.0);
    for (int t = 0; t < n; ++t) {
        if (unit_r_[t] < 0 && unit_d_[t] < 0) throw std::logic_error("subproblem does not cover every taxon");
        w_[t] = 1. / coeff[t];
    }

    x_.assign(static_cast<size_t>(N) * d_, 0.0);
    if (d_ > 0) {
        for (int v = 0; v < N; ++v) {
            double* xv = &x_[static_cast<size_t>(v) * d_];
            const int t = F_.taxon[v];
            if (t >= 0) {
                if (unit_d_[t] >= 0) xv[unit_d_[t]] = w_[t];
            } else {
                for (int e = F_.child_off[v]; e < F_.child_off[v + 1]; ++e) {
                    const double* xc = &x_[static_cast<size_t>(F_.children[e]) * d_];
                    for (int j = 0; j < d_; ++j) xv[j] += xc[j];
                }
            }
        }
    }
    tx_.assign(static_cast<size_t>(F_.k) * d_, 0.0);
    present_.assign(F_.k, {});
    for (int g = 0; g < F_.k; ++g) {
        const int root = F_.tree_off[g + 1] - 1;
        for (int j = 0; j < d_; ++j) {
            tx_[static_cast<size_t>(g) * d_ + j] = x_[static_cast<size_t>(root) * d_ + j];
            if (tx_[static_cast<size_t>(g) * d_ + j] > 0) present_[g].push_back(j);
        }
    }

    a_.assign(N, 0.0); b_.assign(N, 0.0);
    ra_.assign(N, 0); rb_.assign(N, 0);
    ta_.assign(F_.k, 0.0); tb_.assign(F_.k, 0.0);
    tra_.assign(F_.k, 0); trb_.assign(F_.k, 0);
    push_.assign(2 * static_cast<size_t>(N), 0.0);
    dsum_.assign(d_, 0.0);
    xt_.assign(d_, 0.0);
    part_r_.assign(sp.reals.size(), 0);
    part_d_.assign(d_, 0);
}

void Scorer::set_partition(const std::vector<int8_t>& part_r, const std::vector<int8_t>& part_d) {
    part_r_ = part_r;
    part_d_ = part_d;
    const int N = F_.nodes();
    for (int v = 0; v < N; ++v) {                 // children precede parents
        const int t = F_.taxon[v];
        if (t >= 0) {
            if (unit_r_[t] >= 0) {
                const int p = part_r_[unit_r_[t]];
                a_[v] = p == 0; b_[v] = p == 1;
                ra_[v] = p == 0; rb_[v] = p == 1;
            } else {
                const int p = part_d_[unit_d_[t]];
                a_[v] = p == 0 ? w_[t] : 0.0;
                b_[v] = p == 1 ? w_[t] : 0.0;
                ra_[v] = rb_[v] = 0;
            }
        } else {
            double a = 0, b = 0;
            int ra = 0, rb = 0;
            for (int e = F_.child_off[v]; e < F_.child_off[v + 1]; ++e) {
                const int c = F_.children[e];
                a += a_[c]; b += b_[c]; ra += ra_[c]; rb += rb_[c];
            }
            a_[v] = a; b_[v] = b; ra_[v] = ra; rb_[v] = rb;
        }
    }
    for (int g = 0; g < F_.k; ++g) {
        const int root = F_.tree_off[g + 1] - 1;
        ta_[g] = a_[root]; tb_[g] = b_[root];
        tra_[g] = ra_[root]; trb_[g] = rb_[root];
    }
}

double Scorer::score() {
    double s = 0;
    for (int g = 0; g < F_.k; ++g) s += tree_pass(g, false, nullptr, nullptr);
    return s;
}

double Scorer::gains(std::vector<double>& real_gain, std::vector<double>& dummy_gain) {
    real_gain.assign(reals_.size(), 0.0);
    dummy_gain.assign(d_, 0.0);
    double s = 0;
    for (int g = 0; g < F_.k; ++g) s += tree_pass(g, true, real_gain.data(), dummy_gain.data());
    return s;
}

// f at the current generic node (branch values in ga_/gb_/grb_/gx_), with dummy `moved` (or none)
// switched to the other side
double Scorer::generic_value(int D, int moved, PolyEval& pe, const std::vector<int>& pres) {
    int mA = 1, mB = 1;
    for (int j : pres) ((part_d_[j] ^ (j == moved)) == 0 ? mA : mB)++;
    pe.resize(D, mA, mB);
    int ia = 1, ib = 1;
    for (int j : pres) {
        const bool in_a = (part_d_[j] ^ (j == moved)) == 0;
        double* dst = in_a ? &pe.al[static_cast<size_t>(ia++) * D] : &pe.be[static_cast<size_t>(ib++) * D];
        for (int i = 0; i < D; ++i) dst[i] = gx_[i][j];
    }
    for (int i = 0; i < D; ++i) {
        double a = ga_[i], b = gb_[i];
        if (moved >= 0) {
            const double x = gx_[i][moved];
            if (part_d_[moved] == 0) { a -= x; b += x; } else { a += x; b -= x; }
        }
        pe.al[i] = a;
        pe.be[i] = b;
        pe.rho[i] = grb_[i];
    }
    pe.build();
    return pe.eval();
}

// `want` is explicit: with no real (or no dummy) taxa the gain vectors are empty and data() may be null
double Scorer::tree_pass(int g, bool want, double* real_gain, double* dummy_gain) {
    const int off = F_.tree_off[g], end = F_.tree_off[g + 1], root = end - 1;
    const double A = ta_[g], B = tb_[g];
    const int RA = tra_[g], RB = trb_[g];
    const double* X = d_ ? &tx_[static_cast<size_t>(g) * d_] : nullptr;
    const std::vector<int>& pres = present_[g];

    double fsum = 0, G0 = 0, G1 = 0;
    if (want) std::fill(dsum_.begin(), dsum_.end(), 0.0);

    for (int v = off; v < end; ++v) {
        if (F_.taxon[v] >= 0) continue;
        const int cb = F_.child_off[v], ce = F_.child_off[v + 1], nc = ce - cb;
        const double at = A - a_[v], bt = B - b_[v];
        const int rat = RA - ra_[v], rbt = RB - rb_[v];
        const double* xv = d_ ? &x_[static_cast<size_t>(v) * d_] : nullptr;
        for (int j : pres) xt_[j] = X[j] - xv[j];

        if (nc == 2) {
            // branches 0, 1 = children, 2 = top; PA[i] pairs branches (i, i+1), f = sum_i PA[i] PB[i+2]
            const int c0 = F_.children[cb], c1 = F_.children[cb + 1];
            const double a[3] = {a_[c0], a_[c1], at}, b[3] = {b_[c0], b_[c1], bt};
            const int ra3[3] = {ra_[c0], ra_[c1], rat}, rb3[3] = {rb_[c0], rb_[c1], rbt};
            const double* x[3] = {d_ ? &x_[static_cast<size_t>(c0) * d_] : nullptr,
                                  d_ ? &x_[static_cast<size_t>(c1) * d_] : nullptr, xt_.data()};
            double sA[3] = {0, 0, 0}, sB[3] = {double(rb3[0]), double(rb3[1]), double(rb3[2])};
            for (int j : pres) {
                const double x0 = x[0][j], x1 = x[1][j], x2 = x[2][j];
                if (part_d_[j] == 0) { sA[0] += x0 * x1; sA[1] += x1 * x2; sA[2] += x2 * x0; }
                else { sB[0] += x0 * x0; sB[1] += x1 * x1; sB[2] += x2 * x2; }
            }
            double PA[3], PB[3];
            for (int i = 0; i < 3; ++i) {
                PA[i] = a[i] * a[(i + 1) % 3] - sA[i];
                PB[i] = (b[i] * b[i] - sB[i]) / 2;
            }
            const double f = PA[0] * PB[2] + PA[1] * PB[0] + PA[2] * PB[1];
            fsum += f;
            if (!want) continue;

            double dAB[3], dBA[3];
            for (int i = 0; i < 3; ++i) {
                const int j = (i + 1) % 3, l = (i + 2) % 3;
                // the pair opposite branch i is PA[j]; pairs (i,j) and (l,i) lose/gain a_j and a_l
                dAB[i] = ra3[i] > 0 ? PA[j] * b[i] - a[j] * PB[l] - a[l] * PB[j] : 0.0;
                dBA[i] = rb3[i] > 0 ? -PA[j] * (b[i] - 1) + a[j] * PB[l] + a[l] * PB[j] : 0.0;
            }
            push_[2 * c0] = dAB[0] - dAB[2];
            push_[2 * c0 + 1] = dBA[0] - dBA[2];
            push_[2 * c1] = dAB[1] - dAB[2];
            push_[2 * c1 + 1] = dBA[1] - dBA[2];
            G0 += dAB[2];
            G1 += dBA[2];

            for (int jd : pres) {
                const double xs[3] = {x[0][jd], x[1][jd], x[2][jd]};
                const double s = part_d_[jd] == 0 ? -1.0 : 1.0;   // A->B: a -= x, b += x
                double a2[3], b2[3], sA2[3], sB2[3];
                for (int i = 0; i < 3; ++i) {
                    a2[i] = a[i] + s * xs[i];
                    b2[i] = b[i] - s * xs[i];
                    sA2[i] = sA[i] + s * xs[i] * xs[(i + 1) % 3];
                    sB2[i] = sB[i] - s * xs[i] * xs[i];
                }
                double f2 = 0;
                for (int i = 0; i < 3; ++i) {
                    const double pa = a2[i] * a2[(i + 1) % 3] - sA2[i];
                    const double pb = (b2[(i + 2) % 3] * b2[(i + 2) % 3] - sB2[(i + 2) % 3]) / 2;
                    f2 += pa * pb;
                }
                dsum_[jd] += f2 - f;
            }
        } else {
            const int D = nc + 1;
            ga_.resize(D); gb_.resize(D); gra_.resize(D); grb_.resize(D); gx_.resize(D); gd_.resize(2 * D);
            for (int k = 0; k < nc; ++k) {
                const int c = F_.children[cb + k];
                ga_[k] = a_[c]; gb_[k] = b_[c]; gra_[k] = ra_[c]; grb_[k] = rb_[c];
                gx_[k] = d_ ? &x_[static_cast<size_t>(c) * d_] : nullptr;
            }
            ga_[nc] = at; gb_[nc] = bt; gra_[nc] = rat; grb_[nc] = rbt; gx_[nc] = xt_.data();

            const double f0 = generic_value(D, -1, pe_, pres);
            fsum += f0;
            if (!want) continue;

            agg_backup_ = pe_.agg;
            for (int i = 0; i < D; ++i) {
                for (int dir = 0; dir < 2; ++dir) {
                    double delta = 0;
                    if ((dir == 0 ? gra_[i] : grb_[i]) > 0) {
                        const double al0 = pe_.al[i], be0 = pe_.be[i], rho0 = pe_.rho[i];
                        const double s = dir == 0 ? -1.0 : 1.0;   // A->B: a-1, b+1, rho+1
                        pe_.add(i, -1);
                        pe_.al[i] = al0 + s;
                        pe_.be[i] = be0 - s;
                        pe_.rho[i] = rho0 - s;
                        pe_.add(i, 1);
                        delta = pe_.eval() - f0;
                        pe_.al[i] = al0; pe_.be[i] = be0; pe_.rho[i] = rho0;
                        pe_.agg = agg_backup_;
                    }
                    gd_[2 * i + dir] = delta;
                }
            }
            for (int k = 0; k < nc; ++k) {
                const int c = F_.children[cb + k];
                push_[2 * c] = gd_[2 * k] - gd_[2 * nc];
                push_[2 * c + 1] = gd_[2 * k + 1] - gd_[2 * nc + 1];
            }
            G0 += gd_[2 * nc];
            G1 += gd_[2 * nc + 1];
            for (int jd : pres) dsum_[jd] += generic_value(D, jd, pe2_, pres) - f0;
        }
    }

    // per-tree totals: w(PA(g)) w(PB(g)) = all quartets with two taxa on each side (Section 2.5.2)
    double sqA = 0, sqB = 0;
    for (int j : pres) (part_d_[j] == 0 ? sqA : sqB) += X[j] * X[j];
    const double PA = (A * A - RA - sqA) / 2, PB = (B * B - RB - sqB) / 2, tq = PA * PB;
    const double score = 2 * fsum - tq;
    if (!want) return score;

    // real taxa: gain = path sum of pushed deltas + top-branch deltas of every node
    push_[2 * root] = 0;
    push_[2 * root + 1] = 0;
    const double dtot0 = (PA - (A - 1)) * (PB + B) - tq;   // one real taxon A -> B
    const double dtot1 = (PA + A) * (PB - (B - 1)) - tq;   // B -> A
    for (int v = root - 1; v >= off; --v) {
        const int p = F_.parent[v];
        push_[2 * v] += push_[2 * p];
        push_[2 * v + 1] += push_[2 * p + 1];
        const int t = F_.taxon[v];
        if (t >= 0 && unit_r_[t] >= 0) {
            const int i = unit_r_[t];
            real_gain[i] += part_r_[i] == 0 ? 2 * (push_[2 * v] + G0) - dtot0 : 2 * (push_[2 * v + 1] + G1) - dtot1;
        }
    }
    for (int j : pres) {
        const double x = X[j];
        const double dt = part_d_[j] == 0 ? (PA - x * (A - x)) * (PB + x * B) - tq
                                          : (PA + x * A) * (PB - x * (B - x)) - tq;
        dummy_gain[j] += 2 * dsum_[j] - dt;
    }
    return score;
}

void Scorer::swap_real(int i) {
    const int t = reals_[i];
    const int p = part_r_[i];
    for (int g = 0; g < F_.k; ++g) {
        int v = F_.leaf(g, t);
        if (v < 0) continue;
        for (; v >= 0; v = F_.parent[v]) {
            if (p == 0) { a_[v] -= 1; b_[v] += 1; --ra_[v]; ++rb_[v]; }
            else { a_[v] += 1; b_[v] -= 1; ++ra_[v]; --rb_[v]; }
        }
        if (p == 0) { ta_[g] -= 1; tb_[g] += 1; --tra_[g]; ++trb_[g]; }
        else { ta_[g] += 1; tb_[g] -= 1; ++tra_[g]; --trb_[g]; }
    }
    part_r_[i] = static_cast<int8_t>(1 - p);
}

void Scorer::swap_dummy(int j) {
    const int p = part_d_[j];
    const int N = F_.nodes();
    for (int v = 0; v < N; ++v) {
        const double x = x_[static_cast<size_t>(v) * d_ + j];
        if (x == 0) continue;
        if (p == 0) { a_[v] -= x; b_[v] += x; } else { a_[v] += x; b_[v] -= x; }
    }
    for (int g = 0; g < F_.k; ++g) {
        const double x = tx_[static_cast<size_t>(g) * d_ + j];
        if (p == 0) { ta_[g] -= x; tb_[g] += x; } else { ta_[g] += x; tb_[g] -= x; }
    }
    part_d_[j] = static_cast<int8_t>(1 - p);
}

}  // namespace wqfm
