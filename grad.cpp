#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <utility>
#include <vector>
#include <string>
#include <functional>
#include <tr1/unordered_set>
#include <stdexcept>
#include <ctime>
#include <cmath>
#include <fstream>

#ifdef SQLITE_ENABLED
#include <sqlite3.h>
#endif // ifdef SQLITE_ENABLED

using namespace std;

struct Options {
    bool m_stop;
    string m_out;
    double m_convThres, m_damp, m_alphaRel, m_alphaSusp, m_alphaIni, m_alphaMin;

    Options() : m_stop(false), m_convThres(1e-6), m_damp(0.1),
        m_alphaRel(0.95), m_alphaSusp(0.3), m_alphaIni(0.5), m_alphaMin(0.05) { }
};

class Reviewer {
public:
    const string m_name;
    unsigned m_valid;
    double m_alpha[2], m_sqSumEpsilon;

    Reviewer(const char *raw, size_t num, double alphaIni) : m_name(raw, num)
    { m_alpha[0] = m_alpha[1] = alphaIni; }

    double &alpha(int loopCnt) { return m_alpha[loopCnt % 2]; }

    bool operator<(const Reviewer &o) const { return m_name < o.m_name; }

    bool operator==(const Reviewer &o) const { return m_name == o.m_name; }

    void reset() { m_sqSumEpsilon = m_valid = 0; }

    void addEpsilon(double epsilon) {
        m_sqSumEpsilon += epsilon * epsilon;
        m_valid += 1;
    }
};

class PointeeComparator
{
public:
    template<typename T>
    bool operator() (const T *a, const T *b) const { return *a < *b; }
};

class PointeeEqualityPredicate
{
public:
    template<typename T>
    bool operator() (const T *a, const T *b) const { return *a == *b; }
};

static const tr1::hash<string> s_hs;

namespace std {
    namespace tr1 {
        template<>
        struct hash<Reviewer *>
        {
            size_t operator() (const Reviewer *r) const { return s_hs(r->m_name); }
        };
    }
}

// typedef set<Reviewer *, PointeeComparator> People;
typedef tr1::unordered_set<Reviewer *, tr1::hash<Reviewer *>, PointeeEqualityPredicate> People;

double convergence(const People &people)
{
    double n = 0;
    for (People::const_iterator i(people.begin()), iEnd(people.end()); i != iEnd; ++i) {
        const double d = (*i)->alpha(0) - (*i)->alpha(1);
        /* if ( ! isnormal(d)) {
            fprintf(stderr, "%s: %f %f\n", (*i)->m_name.c_str(), (*i)->alpha(0), (*i)->alpha(1));
            if (++errCnt > 10) throw runtime_error("isnormal");
        } */
        n += d * d;
    }
    return n;
}

static const int ERR = -666, NA = -777;

class Restaurant {
public:
    string m_name;
    double m_mu, m_sigma;

    Restaurant(const char *name, unsigned len) : m_name(name, len), m_mu(ERR), m_sigma(ERR) { }
};

namespace std {
    namespace tr1 {
        template<>
        struct hash<Restaurant *>
        {
            size_t operator() (const Restaurant *r) const { return s_hs(r->m_name); }
        };
    }
}

class Review {
public:
    Reviewer *m_reviewer;
    double m_day, m_night;

    static double myAtof(const char * const str)
    {
        /*
        char *cp;
        const double res = strtod(str, &cp);
        return (str == cp) ? NA : res;
        */
        const char *cp = str;
        const int BEFORE = 0;
        double d = 0;
        for (int mode = BEFORE; ; ) {
            const char c = *cp;
            if ('0' <= c && c <= '9') {
                if (mode == BEFORE)
                    d = 10 * d + (c - '0');
                else d += (static_cast<double>(c - '0') / (mode *= 10));
                // printf("%d %.f\n", mode, d);
            } else if (c == '.' && mode == BEFORE) {
                ++mode;
            } else break;
            cp++;
        }
        return (str == cp) ? NA : d;
    }

    static Reviewer *parse(const char * const line, double &day, double &night, double alphaIni) {
        unsigned i = 0;
        Reviewer *reviewer;
        do {
            const char c = line[i];
            if (c == '\0') throw invalid_argument(line);
            if (c == ',') {
                reviewer = new Reviewer(line, i++, alphaIni);
                day = myAtof(line + i);
                break;
            }
        } while (++i);
        do {
            const char c = line[i];
            if (c == '\0') throw invalid_argument(line);
            if (c == ',') {
                night = myAtof(line + i + 1);
                break;
            }
        } while (++i);
        if (day == ERR || night == ERR) throw invalid_argument(line);

        return reviewer;
    }

    Review(const char * const line, People &people, double alphaIni) : m_day(ERR), m_night(ERR) {
        Reviewer * const reviewer = parse(line, m_day, m_night, alphaIni);
        const pair<People::iterator, bool> res(people.insert(reviewer));
        m_reviewer = *res.first;
        if ( ! res.second) {
            delete reviewer;
        } /* else if (people.size() < 30) {
            printf("reviewer inserted: %s\n", m_reviewer->m_reviewer->c_str());
        } */
    }

    double getP() const {
        if (m_day > 0) {
            return (m_night > 0) ? (m_day + m_night) * 0.5 : m_day;
        } else return m_night;
    }
};

typedef vector<Review *> VR;
typedef map<Restaurant *, VR *> Restaurant2Reviews;
// typedef tr1::unordered_map<Restaurant *, VR *> Restaurant2Reviews;

/*
def weighted_incremental_variance(dataWeightPairs):
    sumweight = 0
    mean = 0
    M2 = 0

    for x, weight in dataWeightPairs:  # Alternately "for x, weight in zip(data, weights):"
        temp = weight + sumweight
        delta = x − mean
        R = delta * weight / temp
        mean = mean + R
        if sumweight > 0:
            M2 = M2 + sumweight * delta * R  # Alternatively, "M2 = M2 + weight * delta * (x−mean)"
        sumweight = temp

    variance_n = M2/sumweight
    variance = variance_n * len(dataWeightPairs)/(len(dataWeightPairs) − 1)
*/

void calcMuSigma(const Restaurant2Reviews &rest2rev, int lc)
{
    for (Restaurant2Reviews::const_iterator j(rest2rev.begin()), jEnd(rest2rev.end()); j != jEnd; ++j) {
        const VR &vr = *j->second;
        double mu_j = 0, M2 = 0, sum_alpha_i = 0;
        for (VR::const_iterator ij(vr.begin()), ijEnd(vr.end()); ij != ijEnd; ++ij) {
            const double p_ij = (*ij)->getP(), alpha_i = (*ij)->m_reviewer->alpha(lc);
            const double sum_alpha_i_nx = sum_alpha_i + alpha_i,
                 delta = p_ij - mu_j, R = delta * alpha_i / sum_alpha_i_nx;
            mu_j += R;
            M2 += sum_alpha_i * delta * R;
            sum_alpha_i = sum_alpha_i_nx;
        }
        j->first->m_mu = mu_j;
        const unsigned n = static_cast<unsigned>(vr.size());
        j->first->m_sigma = sqrt(M2 / sum_alpha_i * n / (n - 1));
    }
}

static const double MIN_DELTA = 1e-6;

void calcEpsilon(const Restaurant2Reviews &rest2rev, const People &people)
{
    for (People::const_iterator i(people.begin()), iEnd(people.end()); i != iEnd; ++i) {
        (*i)->reset();
    }

    for (Restaurant2Reviews::const_iterator j(rest2rev.begin()), jEnd(rest2rev.end()); j != jEnd; ++j) {
        const VR &vr = *j->second;
        const double mu_j = j->first->m_mu, sigma_j = j->first->m_sigma;

        for (VR::const_iterator ij(vr.begin()), ijEnd(vr.end()); ij != ijEnd; ++ij) {
            const double p_ij = (*ij)->getP();
            if (p_ij >= 0) {
                const double delta = fabs(p_ij - mu_j),
                     epsilon_ij = (delta < MIN_DELTA) ? 0 : delta / sigma_j;
                (*ij)->m_reviewer->addEpsilon(epsilon_ij);
            }
        }
    }
}

void updateAlpha(const People &people, int lc, double damp, double alphaMin)
{
    // unsigned errCnt = 0;
    for (People::const_iterator i(people.begin()), iEnd(people.end()); i != iEnd; ++i) {
        if ((*i)->m_valid > 0) {
            double alphaNx = exp( - damp * (*i)->m_sqSumEpsilon / (*i)->m_valid);
            /* if ( ! isnormal(alphaNx)) {
                printf("%s: alpha==%f; sum_j epsilon_ij==%f; #epsilon_ij==%u\n",
                    (*i)->m_name.c_str(), (*i)->alpha(lc), (*i)->sqSumEpsilon, (*i)->m_valid);
                if (errCnt++ > 30) throw runtime_error("updateAlpha");
            } */
            (*i)->alpha(lc + 1) = (alphaNx >= alphaMin) ? alphaNx : alphaMin;
        }
    }
}

//output data
void dumpAll(const string &filename, const Restaurant2Reviews &rest2rev, const People &people) {
    ofstream ofs(filename.c_str());

    //regiser shops
    for (Restaurant2Reviews::const_iterator j(rest2rev.begin()), jEnd(rest2rev.end()); j != jEnd; ++j) {
        const Restaurant &r = *(j->first);
        ofs << "shop:" << r.m_name << endl;
    }

    //authorities
    for (People::const_iterator i(people.begin()), iEnd(people.end()); i != iEnd; ++i) {
        const double alpha_i = (*i)->alpha(0);
        ofs << "authority:" << (*i)->m_name << ":" << alpha_i << endl;
    }

    //reviews
    for (Restaurant2Reviews::const_iterator j(rest2rev.begin()), jEnd(rest2rev.end()); j != jEnd; ++j) {
        const Restaurant r = *(j->first);
        const VR &vr = *j->second;

        for (VR::const_iterator ij(vr.begin()), ijEnd(vr.end()); ij != ijEnd; ++ij) {
            const double p_ij = (*ij)->getP();
            if (p_ij >= 0) {
                Reviewer * const i = (*ij)->m_reviewer;
                ofs << "review:" << (r.m_name) << ":" << (i->m_name) << ":" << p_ij << endl;
            }
        }
    }
}

bool calcAuth(const People &people, const Restaurant2Reviews &rest2rev, double criterion, double damp, double minAlha)
{
    int loopCnt = 0, MAX_LOOP = 100;
    for ( ; loopCnt < MAX_LOOP; ++loopCnt) {
        calcMuSigma(rest2rev, loopCnt);
        calcEpsilon(rest2rev, people);
        updateAlpha(people, loopCnt, damp, minAlha);

        const double conv = convergence(people);
        printf("conv == %f\n", conv);
        if (conv < criterion)
            break;
    }

    if (loopCnt >= MAX_LOOP) {
        fprintf(stderr, "Didn't converge after %d times of loop.\n", loopCnt);
        return false;
    } else {
        printf("After %d-th loop:\n", loopCnt);
        return true;
    }
}

void printResult(const People &people, double alphaRel, double alphaSusp)
{
    printf("List of suspicious reviewers are:\n");
    unsigned suspCnt = 0;
    for (People::const_iterator i(people.begin()), iEnd(people.end()); i != iEnd; ++i) {
        const double alpha_i = (*i)->alpha(0);
        if (alpha_i < alphaSusp) {
            printf("%s : %f\n", (*i)->m_name.c_str(), alpha_i);
            ++suspCnt;
        }
    }
    printf("%u suspicious reviewers found.\n", suspCnt);
    printf("On the contrary, ");
    // printf(" list of reliable reviewers are:\n");
    unsigned relCnt = 0;
    for (People::const_iterator i(people.begin()), iEnd(people.end()); i != iEnd; ++i) {
        const double alpha_i = (*i)->alpha(0);
        if (alpha_i > alphaRel) {
            // printf("%s : %f\n", (*i)->m_name.c_str(), alpha_i);
            ++relCnt;
        }
    }
    printf("%u reliable reviewers found.\n", relCnt);
}

void printElapsed(const char *msg, clock_t a, clock_t b)
{
    printf("%s took %.1f sec.\n", msg, double(b - a) / CLOCKS_PER_SEC);
}

int main(int argc, char *argv[]) {
    Options options;

    for (char **arg = argv + 1; *arg; ++arg) {
        if (0 == strcmp("--stop", *arg)) {
            options.m_stop = true;
            ++argv, --argc;
        } else if (0 == strcmp("--out", *arg)) {
            options.m_out = *(++arg);
            argv += 2, argc -= 2;
        } else if (0 == strcmp("--thres", *arg)) {
            options.m_convThres = atof(*(++arg));
            argv += 2, argc -= 2;
        } else if (0 == strcmp("--damp", *arg)) {
            options.m_damp = atof(*(++arg));
            argv += 2, argc -= 2;
        } else if (0 == strcmp("--arel", *arg)) {
            options.m_alphaRel = atof(*(++arg));
            argv += 2, argc -= 2;
        } else if (0 == strcmp("--asusp", *arg)) {
            options.m_alphaSusp = atof(*(++arg));
            argv += 2, argc -= 2;
        } else if (0 == strcmp("--aini", *arg)) {
            options.m_alphaIni = atof(*(++arg));
            argv += 2, argc -= 2;
        } else if (0 == strcmp("--amin", *arg)) {
            options.m_alphaMin = atof(*(++arg));
            argv += 2, argc -= 2;
        }
    }

    if (argc < 2) {
        fprintf(stderr, "USAGE: grad [--stop] [--out output.txt] [--thres 1e-6] [--damp 0.1] input.txt\n");
        exit(10);
    }

#ifdef DEBUG_MYATOF
    for (int i = 1; i < argc; ++i) {
        printf("%s => %f\n", argv[i], Review::myAtof(argv[i]));
    }
#else
    FILE * const fp = fopen(argv[1], "r");
    const size_t buflen = 128;
    char buf[buflen];

    People people;
    Restaurant2Reviews rest2rev;
    unsigned numReviews = 0;
    VR *vr;
    const clock_t t0 = clock();
    while (fgets(buf, buflen, fp)) {
        if (buf[0] == '@') {
            const size_t l = strlen(buf);
            vr = new VR;
            rest2rev.insert(make_pair(new Restaurant(buf + 1, l - 2), vr));
        } else {
            vr->push_back(new Review(buf, people, options.m_alphaIni));
            ++numReviews;
        }
    }
    const clock_t t1 = clock();
    printElapsed("Reading data", t0, t1);
    printf("#people == %u;  #restaurants == %u; numReviews=%u\n",
        static_cast<unsigned>(people.size()), static_cast<unsigned>(rest2rev.size()), numReviews);

    calcAuth(people, rest2rev, options.m_convThres, options.m_damp, options.m_alphaMin);
    const clock_t t2 = clock();
    printElapsed("Authority calculation", t1, t2);

    printResult(people, options.m_alphaRel, options.m_alphaSusp);
    if ( ! options.m_out.empty()) {
        const size_t sz = options.m_out.size();
        if (sz > 4 && string::npos != options.m_out.find(".db", sz - 3)) {
            printf("dumpAllToSqlite()...");
#ifdef SQLITE_ENABLED
            dumpAllToSqlite(options.m_out, rest2rev, people);
#else // ifdef SQLITE_ENABLED
            printf("SQLITE_ENABLED");
            dumpAll(options.m_out, rest2rev, people);
#endif // ifdef SQLITE_ENABLED
        } else
            dumpAll(options.m_out, rest2rev, people);
    }
    const clock_t t3 = clock();
    printElapsed("Output", t2, t3);
#endif // ifdef DEBUG_MYATOF

    if (options.m_stop)
        getchar();

    return 0;
}
