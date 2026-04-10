#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <climits>
#include <cstdio>
#include <cctype>
#include <chrono>

using namespace std;

static int S[5][5];
static constexpr int GAP = 4;
static constexpr int NEG_INF = -1000000000;

static int charIndex(char c) {
    switch (toupper(c)) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return 4;
    }
}

// Sum of all three pairwise scores
static int columnScore(int a, int b, int c) {
    return S[a][b] + S[a][c] + S[b][c];
}

struct FastaEntry {
    string header;
    string seq;
};

static vector<FastaEntry> parseFasta(const string& filename) {
    ifstream fin(filename);
    if (!fin) {
        cerr << "Error: cannot open " << filename << "\n";
        exit(1);
    }

    vector<FastaEntry> entries;
    string line;

    while (getline(fin, line)) {
        if (line.empty()) continue;
        if (line.back() == '\r') line.pop_back();

        if (line[0] == '>') {
            entries.push_back({line.substr(1), ""});
        } else if (!entries.empty()) {
            for (char c : line) {
                if (!isspace((unsigned char)c)) {
                    entries.back().seq += toupper(c);
                }
            }
        }
    }

    return entries;
}

static bool readScoreMatrix(const string& filename) {
    ifstream fin(filename);
    if (!fin) return false;

    string tok;
    for (int j = 0; j < 5; j++) {
        fin >> tok;
    }
    for (int i = 0; i < 5; i++) {
        fin >> tok;
        for (int j = 0; j < 5; j++) {
            fin >> S[i][j];
        }
    }

    return fin.good();
}

// Flat index into 3-D table of size (m+1) x N x K
static inline int idx(int i, int j, int l, int N, int K) {
    return i * N * K + j * K + l;
}

struct Alignment {
    string a, b, c;
    int score;
};

static Alignment align3(const string& seqA, const string& seqB, const string& seqC) {
    const int m = seqA.size();
    const int n = seqB.size();
    const int k = seqC.size();
    const int N = n + 1;
    const int K = k + 1;

    vector<int> dp((m + 1) * N * K, NEG_INF);
    vector<int> trace((m + 1) * N * K, -1);
    dp[0] = 0;

    for (int i = 0; i <= m; i++) {
        int ai = (i > 0) ? charIndex(seqA[i - 1]) : GAP;

        for (int j = 0; j <= n; j++) {
            int bj = (j > 0) ? charIndex(seqB[j - 1]) : GAP;

            for (int l = 0; l <= k; l++) {
                if (i == 0 && j == 0 && l == 0) continue;

                int cl = (l > 0) ? charIndex(seqC[l - 1]) : GAP;
                int best = NEG_INF;
                int dir = -1;

                auto tryMove = [&](int d, int prev, int sc) {
                    int val = prev + sc;
                    if (val > best) {
                        best = val;
                        dir = d;
                    }
                };

                // All 7 transitions — each advances at least one sequence
                if (i > 0 && j > 0 && l > 0) {
                    tryMove(0, dp[idx(i-1, j-1, l-1, N, K)], columnScore(ai, bj, cl));
                }
                if (i > 0 && j > 0) {
                    tryMove(1, dp[idx(i-1, j-1, l, N, K)], columnScore(ai, bj, GAP));
                }
                if (i > 0 && l > 0) {
                    tryMove(2, dp[idx(i-1, j, l-1, N, K)], columnScore(ai, GAP, cl));
                }
                if (j > 0 && l > 0) {
                    tryMove(3, dp[idx(i, j-1, l-1, N, K)], columnScore(GAP, bj, cl));
                }
                if (i > 0) {
                    tryMove(4, dp[idx(i-1, j, l, N, K)], columnScore(ai, GAP, GAP));
                }
                if (j > 0) {
                    tryMove(5, dp[idx(i, j-1, l, N, K)], columnScore(GAP, bj, GAP));
                }
                if (l > 0) {
                    tryMove(6, dp[idx(i, j, l-1, N, K)], columnScore(GAP, GAP, cl));
                }

                dp[idx(i, j, l, N, K)] = best;
                trace[idx(i, j, l, N, K)] = dir;
            }
        }
    }

    // Traceback
    Alignment aln;
    aln.score = dp[idx(m, n, k, N, K)];

    int i = m, j = n, l = k;
    while (i > 0 || j > 0 || l > 0) {
        int dir = trace[idx(i, j, l, N, K)];
        char ca = '-', cb = '-', cc = '-';

        switch (dir) {
            case 0:
                ca = seqA[i - 1]; cb = seqB[j - 1]; cc = seqC[l - 1];
                i--; j--; l--;
                break;
            case 1:
                ca = seqA[i - 1]; cb = seqB[j - 1];
                i--; j--;
                break;
            case 2:
                ca = seqA[i - 1]; cc = seqC[l - 1];
                i--; l--;
                break;
            case 3:
                cb = seqB[j - 1]; cc = seqC[l - 1];
                j--; l--;
                break;
            case 4:
                ca = seqA[i - 1];
                i--;
                break;
            case 5:
                cb = seqB[j - 1];
                j--;
                break;
            case 6:
                cc = seqC[l - 1];
                l--;
                break;
            default:
                i--;
                break;
        }

        aln.a += ca;
        aln.b += cb;
        aln.c += cc;
    }

    reverse(aln.a.begin(), aln.a.end());
    reverse(aln.b.begin(), aln.b.end());
    reverse(aln.c.begin(), aln.c.end());

    return aln;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " input.fasta [score_matrix.txt]\n";
        return 1;
    }

    auto totalStart = chrono::high_resolution_clock::now();

    // Default BLAST-like score matrix
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            if (i == 4 && j == 4)
                S[i][j] = 0;
            else if (i == 4 || j == 4)
                S[i][j] = -8;
            else if (i == j)
                S[i][j] = 5;
            else
                S[i][j] = -4;
        }
    }

    if (argc >= 3 && readScoreMatrix(argv[2])) {
        cout << "Loaded score matrix from " << argv[2] << "\n";
    }

    auto entries = parseFasta(argv[1]);
    if (entries.size() < 3) {
        cerr << "Need at least 3 sequences.\n";
        return 1;
    }

    cout << "Seq 1: " << entries[0].header << "  (len " << entries[0].seq.size() << ")\n";
    cout << "Seq 2: " << entries[1].header << "  (len " << entries[1].seq.size() << ")\n";
    cout << "Seq 3: " << entries[2].header << "  (len " << entries[2].seq.size() << ")\n";

    // Print score matrix
    const char* lab = "ACGT-";
    cout << "\nScore matrix:\n   ";
    for (int j = 0; j < 5; j++) {
        printf(" %3c", lab[j]);
    }
    printf("\n");
    for (int i = 0; i < 5; i++) {
        printf("%c  ", lab[i]);
        for (int j = 0; j < 5; j++) {
            printf(" %3d", S[i][j]);
        }
        printf("\n");
    }

    // Run alignment
    Alignment aln = align3(entries[0].seq, entries[1].seq, entries[2].seq);

    // Compute stats
    int alnLen = aln.a.size();
    int perfectCols = 0;
    for (int i = 0; i < alnLen; i++) {
        if (aln.a[i] != '-' && aln.a[i] == aln.b[i] && aln.b[i] == aln.c[i]) {
            perfectCols++;
        }
    }

    auto totalEnd = chrono::high_resolution_clock::now();
    double totalTime = chrono::duration<double>(totalEnd - totalStart).count();

    cout << "\n================ Results ================\n";
    cout << "Optimal alignment score      : " << aln.score << "\n";
    cout << "Alignment length             : " << alnLen << "\n";
    cout << "Perfectly matched columns    : " << perfectCols << "\n";
    printf("Total time                   : %.3f seconds\n", totalTime);
    cout << "==========================================\n";

    // Write alignment to file
    const char* outFile = "alignment_output.txt";
    ofstream fout(outFile);
    fout << "# Score: " << aln.score << "\n";
    fout << "# Length: " << alnLen << "\n";
    fout << "# Perfect-match columns: " << perfectCols << "\n";

    for (int i = 0; i < 3; i++) {
        fout << ">" << entries[i].header << "\n";

        const string& row = (i == 0) ? aln.a : (i == 1) ? aln.b : aln.c;
        for (int p = 0; p < alnLen; p += 80) {
            fout << row.substr(p, 80) << "\n";
        }
    }

    fout.close();
    cout << "Alignment written to " << outFile << "\n";

    return 0;
}