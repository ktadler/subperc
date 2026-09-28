// subperc_sim: Hamilton backbone plus a repaired pairing, then a D-p scan.
//
// Each instance is a Hamilton cycle on vertices 0..n-1 plus a uniform
// random (d-2)-regular graph on the same vertices. The shortcut graph is a
// configuration model. Self-loops, duplicate pairs, and edges that repeat
// the cycle are removed by random double-edge swaps. The union is a simple
// d-regular graph, written as the cycle plus the shortcuts.
//
// For each D = 0.00, 0.01, ..., 0.79 the program deletes floor(D*n)
// vertices by the floating-range rule: a counter steps by 1/D and the
// vertex floor(counter) is removed. It then scans p. Without --exp it
// records two bracket midpoints, where the largest occupied component
// crosses n^{1/3} (columns p13_*) and where it crosses n^{2/3} (columns
// p23_*). At n = 10^6 the second of these is the cutoff 0.01n used in the
// paper. The figures take the mean of the p23 columns.
//
// Output: <outdir>/table<d>reg.tab
// With --save-instances it also writes <outdir>/graph_<d>_<run>.edges
// ("u v" per line, cycle edges first).
//
// Usage: subperc_sim n d nruns seed outdir [--save-instances] [--exp EXP] [--gap GAP]
// Example: ./subperc_sim 1000000 3 3 12345 out
//   --exp 1/2 or 1/3 records only that cutoff (column p12_* or p13_*).
//   --gap 0.1 does not scan cutoffs. For each D it solves E(D,p)=1±gap and
//   writes the largest component at those two values of p to
//   size<d>_gap<GAP>.tab.
// The run index r = 0,1,... uses the generator seed + 1000003*r + d.
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <unordered_set>

// ---------- union-find with running maximum component size ----------
struct DSU {
    std::vector<int> parent;
    std::vector<int> sz;
    long maxsz = 1;
    void reset(int n) {
        parent.resize(n);
        sz.assign(n, 1);
        for (int i = 0; i < n; i++) parent[i] = i;
        maxsz = 1;
    }
    int find(int x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    }
    void unite(int a, int b) {
        a = find(a); b = find(b);
        if (a == b) return;
        if (sz[a] < sz[b]) std::swap(a, b);
        parent[b] = a;
        sz[a] += sz[b];
        if (sz[a] > maxsz) maxsz = sz[a];
    }
};

static inline long long pair_key(int a, int b, int n) {
    if (a > b) std::swap(a, b);
    return (long long)a * (long long)n + (long long)b;
}

static inline bool ring_adjacent(int u, int v, int n) {
    int diff = u - v;
    if (diff < 0) diff = -diff;
    return diff == 1 || diff == n - 1;
}

// ---------- ring + (d-2)-regular pairing with local repair ----------
// Fills eu/ev with the m2 = n(d-2)/2 shortcut edges. Returns #defects repaired.
static int generate_shortcuts(int n, int d, std::mt19937_64& rng,
                              std::vector<int>& eu, std::vector<int>& ev) {
    const int stars = d - 2;
    const long m2 = (long)n * stars / 2;
    std::vector<int> stubs((size_t)n * stars);
    for (int v = 0; v < n; v++)
        for (int k = 0; k < stars; k++) stubs[(size_t)v * stars + k] = v;
    std::shuffle(stubs.begin(), stubs.end(), rng);

    eu.resize(m2); ev.resize(m2);
    for (long i = 0; i < m2; i++) { eu[i] = stubs[2 * i]; ev[i] = stubs[2 * i + 1]; }

    // classify defects: self-loop, ring collision, duplicate shortcut
    std::unordered_set<long long> exists;
    exists.reserve((size_t)m2 * 2);
    std::vector<long> bad;
    for (long i = 0; i < m2; i++) {
        if (eu[i] == ev[i] || ring_adjacent(eu[i], ev[i], n)) { bad.push_back(i); continue; }
        long long k = pair_key(eu[i], ev[i], n);
        if (exists.count(k)) bad.push_back(i);
        else exists.insert(k);
    }
    std::vector<uint8_t> is_bad(m2, 0);
    for (long i : bad) is_bad[i] = 1;

    std::uniform_int_distribution<long> pick(0, m2 - 1);
    std::uniform_int_distribution<int> coin(0, 1);
    int repaired = 0;
    for (long idx : bad) {
        int u = eu[idx], v = ev[idx];
        for (;;) {
            long j = pick(rng);
            if (j == idx || is_bad[j]) continue;
            int p = eu[j], q = ev[j];
            int nu1, nv1, nu2, nv2;
            if (u == v) {                 // self-loop: (u,u)+(p,q) -> (u,p)+(u,q)
                nu1 = u; nv1 = p; nu2 = u; nv2 = q;
            } else if (coin(rng)) {       // (u,v)+(p,q) -> (u,p)+(v,q)
                nu1 = u; nv1 = p; nu2 = v; nv2 = q;
            } else {                      //            or (u,q)+(v,p)
                nu1 = u; nv1 = q; nu2 = v; nv2 = p;
            }
            if (nu1 == nv1 || nu2 == nv2) continue;
            if (ring_adjacent(nu1, nv1, n) || ring_adjacent(nu2, nv2, n)) continue;
            long long k1 = pair_key(nu1, nv1, n), k2 = pair_key(nu2, nv2, n);
            if (k1 == k2) continue;
            if (exists.count(k1) || exists.count(k2)) continue;
            exists.erase(pair_key(p, q, n));
            exists.insert(k1); exists.insert(k2);
            eu[idx] = nu1; ev[idx] = nv1;
            eu[j] = nu2; ev[j] = nv2;
            is_bad[idx] = 0;
            repaired++;
            break;
        }
    }

    // sanity: every vertex has shortcut-degree d-2, no forbidden edges
    std::vector<int> deg(n, 0);
    for (long i = 0; i < m2; i++) {
        if (eu[i] == ev[i] || ring_adjacent(eu[i], ev[i], n)) {
            fprintf(stderr, "internal error: forbidden edge survived repair\n");
            exit(2);
        }
        deg[eu[i]]++; deg[ev[i]]++;
    }
    for (int v = 0; v < n; v++)
        if (deg[v] != stars) { fprintf(stderr, "internal error: bad degree at %d\n", v); exit(2); }
    return repaired;
}

// ---------- one percolation trial: largest component among alive vertices ----------
static long trial(int n, double p,
                  const std::vector<int>& eu, const std::vector<int>& ev,
                  const std::vector<uint8_t>& deleted,
                  DSU& dsu, std::mt19937_64& rng) {
    dsu.reset(n);
    const bool all = (p >= 1.0);
    const uint64_t thr = all ? UINT64_MAX : (uint64_t)(p * 18446744073709551616.0);
    // ring edges
    for (int u = 0; u < n; u++) {
        int v = (u + 1 == n) ? 0 : u + 1;
        if (deleted[u] || deleted[v]) continue;
        if (all || rng() < thr) dsu.unite(u, v);
    }
    // shortcut edges
    const long m2 = (long)eu.size();
    for (long i = 0; i < m2; i++) {
        if (deleted[eu[i]] || deleted[ev[i]]) continue;
        if (all || rng() < thr) dsu.unite(eu[i], ev[i]);
    }
    return dsu.maxsz;
}

// Bracket-midpoint estimate of the critical p for a given size threshold:
// last fine-grid p with largest < thr, first with largest >= thr; report midpoint.
// start_p is a warm hint from a previous D / coarser criterion.
static double find_pc(int n, double thr, double refine_thr,
                      const std::vector<int>& eu, const std::vector<int>& ev,
                      const std::vector<uint8_t>& deleted,
                      DSU& dsu, std::mt19937_64& rng, double start_p) {
    double p = std::min(std::max(start_p, 0.0), 1.0);
    while (p < 1.0 && (double)trial(n, p, eu, ev, deleted, dsu, rng) < thr)
        p = std::min(p + 0.1, 1.0);
    double step = 0.05;
    for (;;) {
        const long largest = trial(n, p, eu, ev, deleted, dsu, rng);
        if ((double)largest < thr) break;
        if ((double)largest < refine_thr) step = 0.005;
        p -= step;
        if (p <= 0) { p = 0; break; }
    }
    bool saturated = true;
    while (p + 0.005 <= 1.0) {
        if ((double)trial(n, p + 0.005, eu, ev, deleted, dsu, rng) < thr) {
            p += 0.005;
        } else {
            saturated = false;
            break;
        }
    }
    return saturated ? 1.0 : p + 0.0025;
}

// Closed form of the mean offspring. Must match E() in plot_compare.py,
// including the D=0 limit used by p_crit() there.
static double offspring_E(double p, double D, int d) {
    if (p <= 0.0) return 0.0;
    if (p >= 1.0) return 1e300;
    if (D <= 0.0)
        return (d - 3) * p + 2.0 * (d - 2) * p * p / (1.0 - p);
    const int k = (int)std::floor(1.0 / D) - 1;
    const double qk = (k + 2) * D - 1.0;
    const double qk1 = 1.0 - (k + 1) * D;
    const double s = qk * std::pow(p, k) + qk1 * std::pow(p, k + 1);
    return (1.0 - D) * (d - 3) * p
        + 2.0 * p * p * (d - 2) / ((1.0 - p) * (1.0 - p))
            * (1.0 - 2.0 * D - p + p * D + s);
}

// Smallest p in (0,1) with offspring_E >= target, or -1 if the target is
// never reached. E is increasing in p.
static double p_for_E(double D, int d, double target) {
    // Stay a finite distance from p=1: the closed form is 0/0 there and
    // cancellation makes E non-monotone in the last 1e-9.
    const double lo0 = 1e-12, hi0 = 1.0 - 1e-6;
    if (offspring_E(hi0, D, d) < target) return -1.0;
    if (offspring_E(lo0, D, d) >= target) return lo0;
    double lo = lo0, hi = hi0;
    for (int i = 0; i < 80; i++) {
        const double mid = 0.5 * (lo + hi);
        if (offspring_E(mid, D, d) < target) lo = mid;
        else hi = mid;
    }
    return 0.5 * (lo + hi);
}

// Column prefix for a cutoff n^{exp}: p13, p12, p23, or pXX.
static const char* exp_prefix(double e) {
    if (std::fabs(e - 1.0 / 3.0) < 1e-8) return "p13";
    if (std::fabs(e - 0.5) < 1e-8) return "p12";
    if (std::fabs(e - 2.0 / 3.0) < 1e-8) return "p23";
    static char buf[16];
    std::snprintf(buf, sizeof(buf), "p%02d", (int)std::lround(e * 100));
    return buf;
}

int main(int argc, char** argv) {
    if (argc < 6) {
        fprintf(stderr,
                "usage: %s n d nruns seed outdir [--save-instances] [--exp EXP] [--gap GAP]\n",
                argv[0]);
        return 1;
    }
    const int n = atoi(argv[1]);
    const int d = atoi(argv[2]);
    const int nruns = atoi(argv[3]);
    const uint64_t seed = strtoull(argv[4], nullptr, 10);
    const std::string outdir = argv[5];
    bool save_instances = false;
    double exp_only = -1.0;
    double gap = -1.0;
    for (int a = 6; a < argc; a++) {
        if (!strcmp(argv[a], "--save-instances")) save_instances = true;
        else if (!strcmp(argv[a], "--exp") && a + 1 < argc) exp_only = atof(argv[++a]);
        else if (!strcmp(argv[a], "--gap") && a + 1 < argc) gap = atof(argv[++a]);
    }
    if (gap > 0.0 && exp_only > 0.0) {
        fprintf(stderr, "--gap and --exp are different runs\n");
        return 1;
    }
    if (gap > 0.0 && !(gap < 1.0)) {
        fprintf(stderr, "--gap must lie in (0,1)\n");
        return 1;
    }

    if (d < 3) { fprintf(stderr, "need d >= 3\n"); return 1; }
    if (((long)n * (d - 2)) % 2 != 0) { fprintf(stderr, "need n(d-2) even\n"); return 1; }

    const int nD = 80;                       // D = 0.00 .. 0.79
    const double thr13 = std::pow((double)n, 1.0 / 3.0);
    const double thr12 = std::pow((double)n, 0.5);
    const double thr23 = std::pow((double)n, 2.0 / 3.0);
    const double refine_thr = std::pow((double)n, 10.0 / 11.0);
    const bool gap_mode = gap > 0.0;
    const bool single = exp_only > 0.0;
    const double thr_one = single ? std::pow((double)n, exp_only) : 0.0;
    const char* pref_one = single ? exp_prefix(exp_only) : "";
    if (gap_mode)
        fprintf(stderr, "n=%d  gap=%.4f  targets E=%.4f and E=%.4f\n",
                n, gap, 1.0 - gap, 1.0 + gap);
    else if (single)
        fprintf(stderr, "n=%d  cutoff n^{%.4f}=%.0f  columns %s_*\n",
                n, exp_only, thr_one, pref_one);
    else
        fprintf(stderr, "n=%d  n^{1/3}=%.0f  n^{1/2}=%.0f  n^{2/3}=%.0f\n",
                n, thr13, thr12, thr23);

    std::vector<std::vector<double>> p13(nruns, std::vector<double>(nD, 0.0));
    std::vector<std::vector<double>> p12(nruns, std::vector<double>(nD, 0.0));
    std::vector<std::vector<double>> p23(nruns, std::vector<double>(nD, 0.0));
    std::vector<std::vector<double>> pone(nruns, std::vector<double>(nD, 0.0));
    std::vector<double> p_lo(nD, -1.0), p_hi(nD, -1.0);
    std::vector<std::vector<long>> S_lo(nruns, std::vector<long>(nD, -1));
    std::vector<std::vector<long>> S_hi(nruns, std::vector<long>(nD, -1));
    if (gap_mode) {
        for (int ind = 0; ind < nD; ind++) {
            const double D = ind / 100.0;
            p_lo[ind] = p_for_E(D, d, 1.0 - gap);
            p_hi[ind] = p_for_E(D, d, 1.0 + gap);
        }
    }
    std::vector<int> eu, ev;
    std::vector<uint8_t> deleted(n);
    DSU dsu;

    for (int run = 0; run < nruns; run++) {
        std::mt19937_64 rng(seed + 1000003ULL * (uint64_t)run + (uint64_t)d);
        int repaired = generate_shortcuts(n, d, rng, eu, ev);
        fprintf(stderr, "d=%d run=%d: instance ready (%d defects repaired)\n",
                d, run, repaired);

        if (save_instances) {
            std::string fn = outdir + "/graph_" + std::to_string(d) + "_" +
                             std::to_string(run) + ".edges";
            FILE* f = fopen(fn.c_str(), "w");
            if (!f) { fprintf(stderr, "cannot write %s\n", fn.c_str()); return 1; }
            for (int u = 0; u < n; u++) fprintf(f, "%d %d\n", u, (u + 1 == n) ? 0 : u + 1);
            for (size_t i = 0; i < eu.size(); i++) fprintf(f, "%d %d\n", eu[i], ev[i]);
            fclose(f);
        }

        double hint = 0.1;
        for (int ind = 0; ind < nD; ind++) {
            const double D = ind / 100.0;
            std::fill(deleted.begin(), deleted.end(), 0);
            if (D > 0) {                     // floating range: step 1/D along the cycle
                const double k = 1.0 / D;
                for (double counter = 0; counter < n; counter += k)
                    deleted[(int)std::floor(counter)] = 1;
            }
            if (gap_mode) {
                if (p_lo[ind] >= 0.0)
                    S_lo[run][ind] = trial(n, p_lo[ind], eu, ev, deleted, dsu, rng);
                if (p_hi[ind] >= 0.0)
                    S_hi[run][ind] = trial(n, p_hi[ind], eu, ev, deleted, dsu, rng);
                fprintf(stderr, "d=%d run=%d D=%.2f -> S(E=%.2f)=%ld  S(E=%.2f)=%ld\n",
                        d, run, D, 1.0 - gap, S_lo[run][ind], 1.0 + gap, S_hi[run][ind]);
            } else if (single) {
                const double a = find_pc(n, thr_one, refine_thr, eu, ev, deleted, dsu, rng, hint);
                pone[run][ind] = a;
                hint = a;
                fprintf(stderr, "d=%d run=%d D=%.2f -> %s=%.4f\n",
                        d, run, D, pref_one, a);
            } else {
                // n^{1/3}, then n^{1/2}, then n^{2/3}, each warm-started from the last.
                const double a = find_pc(n, thr13, refine_thr, eu, ev, deleted, dsu, rng, hint);
                const double m = find_pc(n, thr12, refine_thr, eu, ev, deleted, dsu, rng, a);
                const double b = find_pc(n, thr23, refine_thr, eu, ev, deleted, dsu, rng, m);
                p13[run][ind] = a;
                p12[run][ind] = m;
                p23[run][ind] = b;
                hint = a;
                fprintf(stderr, "d=%d run=%d D=%.2f -> p13=%.4f  p12=%.4f  p23=%.4f\n",
                        d, run, D, a, m, b);
            }
        }
    }

    char gap_tag[32];
    std::snprintf(gap_tag, sizeof(gap_tag), "%g", gap);
    std::string fn = gap_mode
        ? outdir + "/size" + std::to_string(d) + "_gap" + gap_tag + ".tab"
        : outdir + "/table" + std::to_string(d) + "reg.tab";
    FILE* f = fopen(fn.c_str(), "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", fn.c_str()); return 1; }
    if (gap_mode) {
        fprintf(f, "# n=%d d=%d gap=%g\n", n, d, gap);
        fprintf(f, "D\tp_lo\tp_hi");
        for (int r = 0; r < nruns; r++) fprintf(f, "\tS_lo_%d", r + 1);
        for (int r = 0; r < nruns; r++) fprintf(f, "\tS_hi_%d", r + 1);
        fprintf(f, "\n");
        for (int ind = 0; ind < nD; ind++) {
            fprintf(f, "%.2f", ind / 100.0);
            if (p_lo[ind] < 0.0) fprintf(f, "\tnan");
            else fprintf(f, "\t%.6f", p_lo[ind]);
            if (p_hi[ind] < 0.0) fprintf(f, "\tnan");
            else fprintf(f, "\t%.6f", p_hi[ind]);
            for (int r = 0; r < nruns; r++) {
                if (S_lo[r][ind] < 0) fprintf(f, "\tnan");
                else fprintf(f, "\t%ld", S_lo[r][ind]);
            }
            for (int r = 0; r < nruns; r++) {
                if (S_hi[r][ind] < 0) fprintf(f, "\tnan");
                else fprintf(f, "\t%ld", S_hi[r][ind]);
            }
            fprintf(f, "\n");
        }
        fclose(f);
        fprintf(stderr, "wrote %s\n", fn.c_str());
        return 0;
    }
    fprintf(f, "D");
    if (single) {
        for (int r = 0; r < nruns; r++) fprintf(f, "\t%s_%d", pref_one, r + 1);
    } else {
        for (int r = 0; r < nruns; r++) fprintf(f, "\tp13_%d", r + 1);
        for (int r = 0; r < nruns; r++) fprintf(f, "\tp12_%d", r + 1);
        for (int r = 0; r < nruns; r++) fprintf(f, "\tp23_%d", r + 1);
    }
    fprintf(f, "\n");
    for (int ind = 0; ind < nD; ind++) {
        fprintf(f, "%.2f", ind / 100.0);
        if (single) {
            for (int r = 0; r < nruns; r++) fprintf(f, "\t%.4f", pone[r][ind]);
        } else {
            for (int r = 0; r < nruns; r++) fprintf(f, "\t%.4f", p13[r][ind]);
            for (int r = 0; r < nruns; r++) fprintf(f, "\t%.4f", p12[r][ind]);
            for (int r = 0; r < nruns; r++) fprintf(f, "\t%.4f", p23[r][ind]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
    fprintf(stderr, "wrote %s\n", fn.c_str());
    return 0;
}
