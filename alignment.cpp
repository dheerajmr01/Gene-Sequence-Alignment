/*
 * align3.cpp — Space-Efficient Three-Sequence Global Alignment (Optimized)
 *
 * Algorithm : Hirschberg-style divide-and-conquer on the longest sequence,
 *             keeping two O(m*n) layers (forward / backward) at each level.
 *             Base case (|C| <= 1) uses full 3-D DP with traceback.
 *
 * Optimizations over naive version:
 *   - Flat 1-D arrays for all DP tables (cache-friendly)
 *   - Zero heap allocation in the hot loop (pre-allocated reusable buffers)
 *   - Index ranges instead of substr copies at every recursion level
 *   - Precomputed column-score table (125 entries)
 *   - All mutable state encapsulated in AlignCtx (no globals, reentrant)
 *
 * Space     : O(m*n)   where m,n are lengths of the two shorter sequences
 * Time      : O(m*n*k) where k is the length of the longest sequence
 *
 * Compile   : g++ -O3 -march=native -std=c++17 -o align3 align3.cpp
 * Run       : ./align3 input.fasta [score_matrix.txt]
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <climits>
#include <chrono>
#include <cstdio>
#include <cctype>
#include <cstring>

using namespace std;

// ---------------------------------------------------------------------------
// Score tables
// ---------------------------------------------------------------------------
static int  S[5][5];
static int  CS[5][5][5];      // column score: S[a][b] + S[a][c] + S[b][c]
static constexpr int GAP = 4;
static constexpr int NEG_INF = -1000000000;

static inline int cidx(char c) {
    switch (c) {
        case 'A': case 'a': return 0;  case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;  case 'T': case 't': return 3;
        default: return 4;
    }
}

static void initCS() {
    for (int a = 0; a < 5; a++)
        for (int b = 0; b < 5; b++)
            for (int c = 0; c < 5; c++)
                CS[a][b][c] = S[a][b] + S[a][c] + S[b][c];
}

// ---------------------------------------------------------------------------
// Sequence representation:  integer index array built once, then addressed
// by [start, end) ranges throughout the recursion — no string copies.
// ---------------------------------------------------------------------------
static vector<int> buildIdx(const string& s) {
    vector<int> v(s.size());
    for (size_t i = 0; i < s.size(); i++) v[i] = cidx(s[i]);
    return v;
}

// ---------------------------------------------------------------------------
// FASTA parser
// ---------------------------------------------------------------------------
struct FastaEntry { string header, seq; };

static vector<FastaEntry> parseFasta(const string& fn) {
    ifstream fin(fn);
    if (!fin) { cerr << "Error: cannot open " << fn << "\n"; exit(1); }
    vector<FastaEntry> v;
    string line;
    while (getline(fin, line)) {
        if (line.empty()) continue;
        if (line.back() == '\r') line.pop_back();
        if (line[0] == '>') v.push_back({line.substr(1), ""});
        else if (!v.empty())
            for (char c : line)
                if (!isspace((unsigned char)c)) v.back().seq += toupper(c);
    }
    return v;
}

// ---------------------------------------------------------------------------
// Score-matrix file reader
// ---------------------------------------------------------------------------
static bool readScoreMatrix(const string& fn) {
    ifstream fin(fn);
    if (!fin) return false;
    string tok;
    for (int j = 0; j < 5; j++) fin >> tok;
    for (int i = 0; i < 5; i++) { fin >> tok; for (int j = 0; j < 5; j++) fin >> S[i][j]; }
    return fin.good();
}

// ---------------------------------------------------------------------------
// AlignCtx — encapsulates all mutable alignment state.
// Allocated once per alignment run; passed by reference to all functions.
// Eliminates globals and static locals, making the code reentrant.
// ---------------------------------------------------------------------------
struct AlignCtx {
    // Sequence data
    vector<int> Aidx, Bidx, Cidx;
    string sA, sB, sC;

    // DP layers (two scratch + two result)
    vector<int> layer1, layer2, fwd, bwd;

    // forwardLayer scratch: precomputed per-character score offsets
    // "no C advance" scores (constant across layers, depend only on A/B)
    vector<int> fl_cs_ag, fl_cs_gb, fl_cs_ab;
    // "C advance" scores (rewritten each layer since cl changes)
    vector<int> fl_csC_ag, fl_csC_gb, fl_csC_ab;

    // align3_dc scratch: reversed sub-sequence indices
    vector<int> tmpAr, tmpBr, tmpCr;

    void alloc(int maxM, int maxN, int maxK) {
        int sz = (maxM + 1) * (maxN + 1);
        layer1.resize(sz); layer2.resize(sz);
        fwd.resize(sz);    bwd.resize(sz);

        fl_cs_ag.resize(maxM + 1);
        fl_cs_gb.resize(maxN + 1);
        fl_cs_ab.resize(sz);
        fl_csC_ag.resize(maxM + 1);
        fl_csC_gb.resize(maxN + 1);
        fl_csC_ab.resize(sz);

        tmpAr.resize(maxM);
        tmpBr.resize(maxN);
        tmpCr.resize(maxK);
    }
};

// ---------------------------------------------------------------------------
// Forward-layer computation (flat arrays, zero allocation)
//
// A[0..m-1], B[0..n-1], C[0..k-1]  (index arrays)
// Result written to `out`, sized (m+1)*(n+1).
// Uses `buf1`, `buf2` as scratch (must be >= (m+1)*(n+1)).
// ---------------------------------------------------------------------------
static void forwardLayer(AlignCtx& ctx,
                         const int* A, int m,
                         const int* B, int n,
                         const int* C, int k,
                         int* out,
                         int* buf1,
                         int* buf2)
{
    const int W = n + 1;           // row width
    int* prev = buf1;
    int* curr = buf2;

    /* layer l = 0: pairwise A vs B, C = all gaps */
    prev[0] = 0;
    for (int j = 1; j <= n; j++)
        prev[j] = prev[j-1] + CS[GAP][B[j-1]][GAP];
    for (int i = 1; i <= m; i++) {
        prev[i * W] = prev[(i-1) * W] + CS[A[i-1]][GAP][GAP];
        const int ai = A[i - 1];
        for (int j = 1; j <= n; j++) {
            const int bj = B[j - 1];
            int v = prev[(i-1)*W + j-1] + CS[ai][bj][GAP];    // diag
            int w;
            w = prev[(i-1)*W + j] + CS[ai][GAP][GAP]; if (w > v) v = w;  // up
            w = prev[i*W + j-1]   + CS[GAP][bj][GAP];  if (w > v) v = w;  // left
            prev[i*W + j] = v;
        }
    }

    /* Precompute per-character score offsets for the "not advancing C" moves.
       These depend only on A[i]/B[j] and are constant across all layers.
       Buffers are preallocated in ctx — no allocation here. */
    int* cs_ag = ctx.fl_cs_ag.data();
    int* cs_gb = ctx.fl_cs_gb.data();
    int* cs_ab = ctx.fl_cs_ab.data();
    for (int i = 1; i <= m; i++) cs_ag[i] = CS[A[i-1]][GAP][GAP];
    for (int j = 1; j <= n; j++) cs_gb[j] = CS[GAP][B[j-1]][GAP];
    for (int i = 1; i <= m; i++)
        for (int j = 1; j <= n; j++)
            cs_ab[i*W+j] = CS[A[i-1]][B[j-1]][GAP];

    /* Pointers into ctx buffers for "advancing C" scores (rewritten each layer) */
    int* csC_ag = ctx.fl_csC_ag.data();
    int* csC_gb = ctx.fl_csC_gb.data();
    int* csC_ab = ctx.fl_csC_ab.data();

    /* layers l = 1 … k */
    for (int l = 1; l <= k; l++) {
        const int cl = C[l - 1];

        /* Recompute per-character scores for "advancing C" moves this layer.
           Only the values change (cl differs each iteration); no resize needed. */
        const int csC_gg = CS[GAP][GAP][cl];
        for (int i = 1; i <= m; i++) csC_ag[i] = CS[A[i-1]][GAP][cl];
        for (int j = 1; j <= n; j++) csC_gb[j] = CS[GAP][B[j-1]][cl];
        for (int i = 1; i <= m; i++)
            for (int j = 1; j <= n; j++)
                csC_ab[i*W+j] = CS[A[i-1]][B[j-1]][cl];

        /* i=0, j=0 */
        curr[0] = prev[0] + csC_gg;

        /* i=0, j>0  (only 3 transitions: prev->advance C with gap/B, curr->B only) */
        for (int j = 1; j <= n; j++) {
            int v = prev[j] + csC_gg;                              // gap,gap,C
            int w;
            w = prev[j-1] + csC_gb[j];       if (w > v) v = w;    // gap,B,C
            w = curr[j-1] + cs_gb[j];        if (w > v) v = w;    // gap,B,gap
            curr[j] = v;
        }

        /* i>0 rows */
        for (int i = 1; i <= m; i++) {
            const int row  = i * W;
            const int prow = (i-1) * W;

            /* j=0 column (only 3 transitions) */
            {
                int v = prev[row] + csC_gg;                                // gap,gap,C
                int w;
                w = prev[prow] + csC_ag[i];         if (w > v) v = w;     // A,gap,C
                w = curr[prow] + cs_ag[i];           if (w > v) v = w;     // A,gap,gap
                curr[row] = v;
            }

            /* j>0: all 7 transitions valid */
            const int ag_noC = cs_ag[i];
            const int ag_C   = csC_ag[i];
            for (int j = 1; j <= n; j++) {
                int v, w;

                /* 4 transitions advancing C (from prev) */
                v = prev[prow + j-1] + csC_ab[row+j];     // A,B,C
                w = prev[prow + j]   + ag_C;               if (w > v) v = w;  // A,gap,C
                w = prev[row  + j-1] + csC_gb[j];          if (w > v) v = w;  // gap,B,C
                w = prev[row  + j]   + csC_gg;             if (w > v) v = w;  // gap,gap,C

                /* 3 transitions NOT advancing C (from curr) */
                w = curr[prow + j-1] + cs_ab[row+j];       if (w > v) v = w;  // A,B,gap
                w = curr[prow + j]   + ag_noC;              if (w > v) v = w;  // A,gap,gap
                w = curr[row  + j-1] + cs_gb[j];           if (w > v) v = w;  // gap,B,gap

                curr[row + j] = v;
            }
        }
        swap(prev, curr);
    }
    memcpy(out, prev, (size_t)(m+1) * W * sizeof(int));
}

// ---------------------------------------------------------------------------
// 3-D DP with traceback (base case, |C| <= 1)
// ---------------------------------------------------------------------------
struct Alignment { string a, b, c; int score = 0; };

static Alignment dp3d(const AlignCtx& ctx,
                      const int* A, int m,
                      const int* B, int n,
                      const int* C, int k,
                      int aOff, int bOff, int cOff)
{
    const int K1 = k + 1, W = (n+1)*K1;
    vector<int>         dp((m+1)*W, NEG_INF);
    vector<signed char> tr((m+1)*W, -1);

    dp[0] = 0;

    for (int i = 0; i <= m; i++) {
        for (int j = 0; j <= n; j++) {
            for (int l = 0; l <= k; l++) {
                if (i==0 && j==0 && l==0) continue;
                const int idx = i*W + j*K1 + l;
                int best = NEG_INF; signed char dir = -1; int v;

                #define TRY(d,x) { v=(x); if(v>best){best=v;dir=(d);} }
                if (i>0&&j>0&&l>0) TRY(0, dp[(i-1)*W+(j-1)*K1+l-1] + CS[A[i-1]][B[j-1]][C[l-1]])
                if (i>0&&j>0)      TRY(1, dp[(i-1)*W+(j-1)*K1+l]   + CS[A[i-1]][B[j-1]][GAP])
                if (i>0&&l>0)      TRY(2, dp[(i-1)*W+j*K1+l-1]     + CS[A[i-1]][GAP][C[l-1]])
                if (j>0&&l>0)      TRY(3, dp[i*W+(j-1)*K1+l-1]     + CS[GAP][B[j-1]][C[l-1]])
                if (i>0)           TRY(4, dp[(i-1)*W+j*K1+l]       + CS[A[i-1]][GAP][GAP])
                if (j>0)           TRY(5, dp[i*W+(j-1)*K1+l]       + CS[GAP][B[j-1]][GAP])
                if (l>0)           TRY(6, dp[i*W+j*K1+l-1]         + CS[GAP][GAP][C[l-1]])
                #undef TRY

                dp[idx] = best;
                tr[idx] = dir;
            }
        }
    }

    Alignment aln;
    aln.score = dp[m*W + n*K1 + k];
    int i=m, j=n, l=k;
    while (i>0 || j>0 || l>0) {
        char ca='-', cb='-', cc='-';
        switch (tr[i*W+j*K1+l]) {
            case 0: ca=ctx.sA[aOff+i-1];cb=ctx.sB[bOff+j-1];cc=ctx.sC[cOff+l-1];--i;--j;--l;break;
            case 1: ca=ctx.sA[aOff+i-1];cb=ctx.sB[bOff+j-1];--i;--j;break;
            case 2: ca=ctx.sA[aOff+i-1];cc=ctx.sC[cOff+l-1];--i;--l;break;
            case 3: cb=ctx.sB[bOff+j-1];cc=ctx.sC[cOff+l-1];--j;--l;break;
            case 4: ca=ctx.sA[aOff+i-1];--i;break;
            case 5: cb=ctx.sB[bOff+j-1];--j;break;
            case 6: cc=ctx.sC[cOff+l-1];--l;break;
            default: --i; break;
        }
        aln.a += ca; aln.b += cb; aln.c += cc;
    }
    reverse(aln.a.begin(), aln.a.end());
    reverse(aln.b.begin(), aln.b.end());
    reverse(aln.c.begin(), aln.c.end());
    return aln;
}

// ---------------------------------------------------------------------------
// Divide-and-conquer  (no string copies — index ranges only)
// ---------------------------------------------------------------------------
static void align3_dc(AlignCtx& ctx,
                      int aS, int aE,    // A range [aS, aE)
                      int bS, int bE,    // B range
                      int cS, int cE,    // C range (divide dimension)
                      string& alnA, string& alnB, string& alnC)
{
    const int m = aE - aS, n = bE - bS, k = cE - cS;
    const int W = n + 1;

    /* ---------- base case ---------- */
    if (k <= 1) {
        Alignment a = dp3d(ctx,
                           ctx.Aidx.data()+aS, m,
                           ctx.Bidx.data()+bS, n,
                           ctx.Cidx.data()+cS, k,
                           aS, bS, cS);
        alnA = move(a.a);  alnB = move(a.b);  alnC = move(a.c);
        return;
    }

    /* ---------- divide ---------- */
    const int mid = k / 2;   // C_left = [cS, cS+mid),  C_right = [cS+mid, cE)
    const int kR  = k - mid;

    /* forward over C_left */
    forwardLayer(ctx,
                 ctx.Aidx.data()+aS, m,
                 ctx.Bidx.data()+bS, n,
                 ctx.Cidx.data()+cS, mid,
                 ctx.fwd.data(), ctx.layer1.data(), ctx.layer2.data());

    /* backward over C_right:
       The backward pass reverses A, B, and C_right, then runs forwardLayer.
       This is correct because CS[a][b][c] = S[a][b] + S[a][c] + S[b][c]
       depends only on character identities, not positions. The 7-transition
       recurrence is symmetric under simultaneous reversal of all three
       sequences: each transition (e.g. advance A+B+C) has an identical mirror.
       Therefore: optimal_score = max_{i,j} Fwd(i,j,mid) + Fwd_rev(m-i,n-j,kR)
       where Fwd_rev runs forward on the reversed right-half sequences. */
    {
        int* tmpAr = ctx.tmpAr.data();
        int* tmpBr = ctx.tmpBr.data();
        int* tmpCr = ctx.tmpCr.data();
        for (int i = 0; i < m; i++)  tmpAr[i] = ctx.Aidx[aE - 1 - i];
        for (int i = 0; i < n; i++)  tmpBr[i] = ctx.Bidx[bE - 1 - i];
        for (int i = 0; i < kR; i++) tmpCr[i] = ctx.Cidx[cE - 1 - i];

        forwardLayer(ctx,
                     tmpAr, m,
                     tmpBr, n,
                     tmpCr, kR,
                     ctx.bwd.data(), ctx.layer1.data(), ctx.layer2.data());
    }

    /* find optimal split (i*, j*) */
    int bestI = 0, bestJ = 0;
    long long bestVal = LLONG_MIN;
    for (int i = 0; i <= m; i++) {
        for (int j = 0; j <= n; j++) {
            long long v = (long long)ctx.fwd[i*W + j]
                        + (long long)ctx.bwd[(m-i)*W + (n-j)];
            if (v > bestVal) { bestVal = v; bestI = i; bestJ = j; }
        }
    }

    /* ---------- conquer ---------- */
    string la, lb, lc, ra, rb, rc;
    align3_dc(ctx, aS, aS+bestI,  bS, bS+bestJ,  cS, cS+mid,   la, lb, lc);
    align3_dc(ctx, aS+bestI, aE,  bS+bestJ, bE,  cS+mid, cE,   ra, rb, rc);

    alnA = std::move(la);  alnA.append(ra);
    alnB = std::move(lb);  alnB.append(rb);
    alnC = std::move(lc);  alnC.append(rc);
}

// ---------------------------------------------------------------------------
// Peak memory (Linux)
// ---------------------------------------------------------------------------
static double getPeakMemMB() {
    ifstream f("/proc/self/status");
    string line;
    while (getline(f, line))
        if (line.rfind("VmPeak:", 0) == 0) {
            long long kb = 0; sscanf(line.c_str(), "VmPeak: %lld kB", &kb);
            return kb / 1024.0;
        }
    return -1;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " input.fasta [score_matrix.txt]\n\n"
             << "Score-matrix file format (5x5, header row):\n"
             << "    A   C   G   T   -\n"
             << "A   5  -4  -4  -4  -8\n"
             << "C  -4   5  -4  -4  -8\n"
             << "G  -4  -4   5  -4  -8\n"
             << "T  -4  -4  -4   5  -8\n"
             << "-  -8  -8  -8  -8   0\n";
        return 1;
    }

    /* ---- default BLAST scores ---- */
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++) {
            if (i==4 && j==4) S[i][j]=0;
            else if (i==4 || j==4) S[i][j]=-8;
            else if (i==j) S[i][j]=5;
            else S[i][j]=-4;
        }

    if (argc >= 3) {
        if (!readScoreMatrix(argv[2]))
            cerr << "Warning: could not read " << argv[2] << "; using defaults.\n";
        else
            cout << "Loaded score matrix from " << argv[2] << "\n";
    }
    initCS();

    /* ---- read sequences ---- */
    auto entries = parseFasta(argv[1]);
    if (entries.size() < 3) { cerr << "Need >= 3 sequences.\n"; return 1; }

    cout << "Seq 1: " << entries[0].header << "  (len " << entries[0].seq.size() << ")\n";
    cout << "Seq 2: " << entries[1].header << "  (len " << entries[1].seq.size() << ")\n";
    cout << "Seq 3: " << entries[2].header << "  (len " << entries[2].seq.size() << ")\n";

    /* ---- assign C = longest ---- */
    int order[3] = {0, 1, 2};
    string seqs[3] = {entries[0].seq, entries[1].seq, entries[2].seq};
    {
        int mx = 2;
        if (seqs[0].size() >= seqs[1].size() && seqs[0].size() >= seqs[2].size()) mx = 0;
        else if (seqs[1].size() >= seqs[0].size() && seqs[1].size() >= seqs[2].size()) mx = 1;
        if (mx != 2) { swap(seqs[mx], seqs[2]); swap(order[mx], order[2]); }
        // Also ensure seqs[0].size() <= seqs[1].size() for slightly better cache behavior
        if (seqs[0].size() > seqs[1].size()) { swap(seqs[0], seqs[1]); swap(order[0], order[1]); }
    }

    /* ---- initialize alignment context ---- */
    AlignCtx ctx;
    ctx.sA = seqs[0]; ctx.sB = seqs[1]; ctx.sC = seqs[2];
    ctx.Aidx = buildIdx(ctx.sA);
    ctx.Bidx = buildIdx(ctx.sB);
    ctx.Cidx = buildIdx(ctx.sC);

    const int m = (int)ctx.sA.size(), n = (int)ctx.sB.size(), k = (int)ctx.sC.size();

    cout << "\nD&C splits sequence of length " << k
         << "\nLayer size: " << m << " x " << n
         << " = " << (long long)(m+1)*(n+1) << " cells\n";

    /* print score matrix */
    const char* lab = "ACGT-";
    cout << "\nScore matrix:\n   ";
    for (int j = 0; j < 5; j++) printf(" %3c", lab[j]);  printf("\n");
    for (int i = 0; i < 5; i++) { printf("%c  ", lab[i]); for (int j = 0; j < 5; j++) printf(" %3d", S[i][j]); printf("\n"); }

    /* ---- allocate reusable buffers ---- */
    ctx.alloc(m, n, k);

    /* ---- run ---- */
    auto t1 = chrono::high_resolution_clock::now();

    string alnA, alnB, alnC;
    align3_dc(ctx, 0, m, 0, n, 0, k, alnA, alnB, alnC);

    auto t2 = chrono::high_resolution_clock::now();
    double elapsed = chrono::duration<double>(t2 - t1).count();

    /* ---- reorder to original input order ---- */
    string outAln[3];
    outAln[order[0]] = alnA;
    outAln[order[1]] = alnB;
    outAln[order[2]] = alnC;

    /* ---- statistics ---- */
    const int alnLen = (int)outAln[0].size();
    int score = 0, perfectCols = 0;
    for (int i = 0; i < alnLen; i++) {
        score += CS[cidx(outAln[0][i])][cidx(outAln[1][i])][cidx(outAln[2][i])];
        if (outAln[0][i] != '-' && outAln[0][i] == outAln[1][i] && outAln[1][i] == outAln[2][i])
            perfectCols++;
    }

    double peakMB = getPeakMemMB();

    cout << "\n================ Results ================\n";
    cout << "Optimal alignment score      : " << score   << "\n";
    cout << "Alignment length             : " << alnLen  << "\n";
    cout << "Perfectly matched columns    : " << perfectCols << "\n";
    printf("Running time                 : %.3f seconds\n", elapsed);
    if (peakMB > 0) printf("Peak memory (VmPeak)         : %.1f MB\n", peakMB);
    cout << "==========================================\n";

    /* ---- conserved regions ---- */
    cout << "\n--- Conserved Region Analysis ---\n";
    constexpr int WINDOW = 10;
    vector<pair<int,int>> regions;
    int runStart = -1;
    for (int i = 0; i < alnLen; i++) {
        bool pm = outAln[0][i] != '-' && outAln[0][i] == outAln[1][i] && outAln[1][i] == outAln[2][i];
        if (pm) { if (runStart < 0) runStart = i; }
        else { if (runStart >= 0 && i - runStart >= WINDOW) regions.push_back({runStart, i-1}); runStart = -1; }
    }
    if (runStart >= 0 && alnLen - runStart >= WINDOW) regions.push_back({runStart, alnLen-1});

    if (regions.empty()) cout << "No conserved blocks of >= " << WINDOW << " consecutive perfect-match columns.\n";
    else {
        cout << regions.size() << " conserved block(s) of >= " << WINDOW << " perfect-match columns:\n";
        for (auto& [s,e] : regions) printf("  columns %5d - %5d   (length %d)\n", s+1, e+1, e-s+1);
    }

    /* ---- write alignment ---- */
    const char* outFile = "alignment_output.txt";
    ofstream fout(outFile);
    fout << "# Score: " << score << "\n# Length: " << alnLen
         << "\n# Perfect-match columns: " << perfectCols << "\n";
    for (int i = 0; i < 3; i++) {
        fout << ">" << entries[i].header << "\n";
        for (int p = 0; p < alnLen; p += 80) fout << outAln[i].substr(p, 80) << "\n";
    }
    fout.close();
    cout << "\nFull alignment written to " << outFile << "\n";

    return 0;
}
