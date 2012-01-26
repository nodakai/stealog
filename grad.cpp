#include <cstdio>
#include <cstdlib>
#include <map>
#include <utility>
#include <vector>
#include <string>
#include <functional>
#include <set>
#include <stdexcept>
#include <ctime>
#include <cmath>

using namespace std;

class Reviewer {
public:
    static const double ALPHA_INIT;

    const string m_name;
    unsigned m_valid;
    double m_alpha[2], sqSumEpsilon;

    Reviewer(const char *raw, size_t num) : m_name(raw, num)
    { m_alpha[0] = m_alpha[1] = ALPHA_INIT; }

    double &alpha(int loopCnt) { return m_alpha[loopCnt % 2]; }

    bool operator<(const Reviewer &o) const {
        return m_name < o.m_name;
    }
};

const double Reviewer::ALPHA_INIT = 0.5;

typedef set<Reviewer *> People;

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
    double mu, sigma;

    Restaurant(const char *name) : m_name(name), mu(ERR), sigma(ERR) { }
};

class Review {
public:
    Reviewer *m_reviewer;
    double m_day, m_night;

    static double myAtof(const char * const str)
    {
        const char *cp = str;
        int BEFORE = -1, mode = BEFORE;
        double d = 0;
        while (true) {
            const char c = *cp;
            if ('0' <= c && c <= '9') {
                d = 10 * d + (c - '0');
                if (mode > BEFORE)
                    ++mode;
            } else if (c == '.' && mode == BEFORE) {
                ++mode;
            } else break;
            cp++;
        }
        while (mode-- > 0) {
            d /= 10;
        }
        return (str == cp) ? NA : d;
    }

    Reviewer *parse(const char * const line, double &day, double &night) {
        int i = 0;
        Reviewer *reviewer;
        do {
            const char c = line[i];
            if (c == '\0') throw invalid_argument(line);
            if (c == ',') {
                reviewer = new Reviewer(line, i++);
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

    Review(const char * const line, People &people) : m_day(ERR), m_night(ERR) {
        Reviewer * const reviewer = parse(line, m_day, m_night);
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
        j->first->mu = mu_j;
        const int n = vr.size();
        j->first->sigma = sqrt(M2 / sum_alpha_i * n / (n - 1));
    }
}

static const double MIN_DELTA = 1e-6;

void calcEpsilon(const Restaurant2Reviews &rest2rev, const People &people)
{
    for (People::const_iterator i(people.begin()), iEnd(people.end()); i != iEnd; ++i) {
        (*i)->sqSumEpsilon = (*i)->m_valid = 0;
    }

    for (Restaurant2Reviews::const_iterator j(rest2rev.begin()), jEnd(rest2rev.end()); j != jEnd; ++j) {
        const VR &vr = *j->second;
        const double mu_j = j->first->mu, sigma_j = j->first->sigma;

        for (VR::const_iterator ij(vr.begin()), ijEnd(vr.end()); ij != ijEnd; ++ij) {
            const double p_ij = (*ij)->getP();
            if (p_ij >= 0) {
                const double delta = fabs(p_ij - mu_j),
                     epsilon_ij = (delta < MIN_DELTA) ? 0 : delta / sigma_j;
                Reviewer * const i = (*ij)->m_reviewer;
                i->sqSumEpsilon -= epsilon_ij * epsilon_ij;
                i->m_valid += 1;
            }
        }
    }
}

void updateAlpha(const People &people, int lc)
{
    unsigned errCnt = 0;
    for (People::const_iterator i(people.begin()), iEnd(people.end()); i != iEnd; ++i) {
        if ((*i)->m_valid > 0) {
            double alphaNx = exp((*i)->sqSumEpsilon / (*i)->m_valid * 0.1);
            if ( ! isnormal(alphaNx)) {
                printf("%s: alpha==%f; sum_j epsilon_ij==%f; #epsilon_ij==%u\n", (*i)->m_name.c_str(), (*i)->alpha(lc), (*i)->sqSumEpsilon, (*i)->m_valid);
                if (errCnt++ > 30) throw runtime_error("updateAlpha");
            }
            if ( ! (alphaNx > 0.05))
                alphaNx = 0.05;
            (*i)->alpha(lc + 1) = alphaNx;
        }
    }
}

void calcAuth(const People &people, const Restaurant2Reviews &rest2rev)
{
    int loopCnt = 0, MAX_LOOP = 100;
    for ( ; loopCnt < MAX_LOOP; ++loopCnt) {
        calcMuSigma(rest2rev, loopCnt);
        calcEpsilon(rest2rev, people);
        updateAlpha(people, loopCnt);

        const double conv = convergence(people);
        printf("conv == %f\n", conv);
        if (conv < 0.05)
            break;
    }

    if (loopCnt >= MAX_LOOP) {
        fprintf(stderr, "Didn't converge after %d times of loop.\n", loopCnt);
    } else {
        printf("After %d-th loop:\n", loopCnt);
        printf("List of suspicious reviewers are:\n");
        unsigned suspCnt = 0;
        for (People::const_iterator i(people.begin()), iEnd(people.end()); i != iEnd; ++i) {
            const double alpha_i = (*i)->alpha(loopCnt);
            if (alpha_i < 0.3) {
                printf("%s : %f\n", (*i)->m_name.c_str(), alpha_i);
                ++suspCnt;
            }
        }
        printf("%u suspicious reviewers found.\n", suspCnt);
    }
}

void printElapsed(const char *msg, clock_t a, clock_t b)
{
    printf("%s took %.1f sec.\n", msg, double(b - a) / CLOCKS_PER_SEC);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "USAGE: grad input.txt\n");
        exit(10);
    }

    /* for (int i = 1; i < argc || exit(10); ++i) {
        printf("%s => %f\n", argv[i], Review::myAtof(argv[i]));
    }
    return 0; */

    FILE * const fp = fopen(argv[1], "r");
    const size_t buflen = 128;
    char buf[buflen];

    People people;
    Restaurant2Reviews rest2rev;
    unsigned count = 0;
    VR *vr;
    const clock_t t0 = clock();
    while (fgets(buf, buflen, fp)) {
        if (buf[0] == '@') {
            Restaurant *rest = new Restaurant(buf + 1);
            vr = new VR;
            rest2rev.insert(make_pair(rest, vr));
        } else {
            vr->push_back(new Review(buf, people));
            ++count;
        }
    }
    const clock_t t1 = clock();
    printElapsed("Reading data", t0, t1);
    printf("#people == %u;  #rest2rev == %u; count=%u\n",
        unsigned(people.size()), unsigned(rest2rev.size()), count);

    calcAuth(people, rest2rev);
    const clock_t t2 = clock();
    printElapsed("Authority calculation", t1, t2);

    return 0;
}
