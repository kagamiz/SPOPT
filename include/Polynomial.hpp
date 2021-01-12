#ifndef __POLYNOMIAL_HPP__
#define __POLYNOMIAL_HPP__

#include "StandardLibraries.hpp"

namespace SPOPT {
    using Term = std::vector<int>;

    class Monomial {
        public:
            Monomial() : coefficient(0), degree(0) {}
            Monomial(Term _t, double _coef, bool sorted = false) : term(_t), coefficient(_coef), degree(_t.size()) { if (!sorted) std::sort(term.begin(), term.end()); }
            Monomial(double _coef) : term({}), coefficient(_coef), degree(0) {}
            Monomial(Term _t, bool sorted = false) : term(_t), coefficient(1), degree(_t.size()) { if (!sorted) std::sort(term.begin(), term.end()); }

            Monomial operator * (const Monomial &rhs) {
                Term newTerm;
                std::merge(this->term.begin(), this->term.end(), rhs.term.begin(), rhs.term.end(), std::back_inserter(newTerm));
                return Monomial(newTerm, this->coefficient * rhs.coefficient, /* sorted = */true);
            }

            Monomial &operator *= (const Monomial &rhs) {
                *this = *this * rhs;
                return *this;
            }

            Monomial operator * (const double &c) {
                return *this * Monomial(c);
            }

            Monomial &operator *= (const double &c) {
                *this *= Monomial(c);
                return *this;
            }

            Monomial operator * (const Term &t) {
                return *this * Monomial(t);
            }

            Monomial &operator *= (const Term &t) {
                *this *= Monomial(t);
                return *this;
            }

            Monomial DifferentiateBy(int ind) {
                int ct = std::count_if(this->term.begin(), this->term.end(), [&](int x) { return x == ind; });    
                if (ct == 0) {
                    return Monomial();
                }
                Term t;
                bool skipped = false;
                for (auto var : this->term) {
                    if (var != ind || skipped) {
                        t.emplace_back(var);
                    }
                    else {
                        skipped = true;
                    }
                }
                return Monomial(t, this->coefficient * ct, /* sorted = */true);
            }

            std::string ToString(bool oneIndexed = false, bool withSign = true) {
                std::ostringstream oss;
                if (withSign) {
                    oss << std::showpos << coefficient;
                    oss << std::noshowpos;
                }
                else {
                    oss << coefficient;
                }
                for (int i = 0; i < term.size();) {
                    int cnt = 1;
                    int j;
                    for (j = i + 1; j < term.size() && term[i] == term[j]; j++) {
                        cnt++;
                    }
                    oss << "*x[" << term[i] + oneIndexed << "]";
                    if (cnt >= 2) {
                        oss << "^" << cnt;
                    }
                    i = j; 
                }
                return oss.str();
            }

            double Evaluate(std::vector<double> &values)
            {
                double ret = coefficient;
                for (auto &ind : term) {
                    if (ind >= values.size()) return 0;
                    ret *= values[ind];
                }
                return ret;
            }

            int degree;
            Term term;
            double coefficient;
    };

    class Polynomial {
        public:
            Polynomial()
            { 
                degree = 0; maxIndex = -1;
            }

            Polynomial(Monomial m)
            {
                monomials[m.term] = m.coefficient;
                degree = m.degree;
                maxIndex = -1;
                if (m.term.size()) {
                    maxIndex = *max_element(m.term.begin(), m.term.end());
                }
            }

            Polynomial(std::map<Term, double> &m) : monomials(m)
            {
                degree = 0;
                maxIndex = -1;
                for (auto monomial : monomials) {
                    degree = std::max(degree, (int)monomial.first.size());
                    if (monomial.first.size()) {
                        maxIndex = *max_element(monomial.first.begin(), monomial.first.end());
                    }
                }
            }

            Polynomial(std::vector<Monomial> &m)
            {
                degree = 0;
                maxIndex = -1;
                for (auto monomial : m) {
                    if (monomials.find(monomial.term) == monomials.end()) {
                        monomials[monomial.term] = monomial.coefficient;
                    }
                    else {
                        monomials[monomial.term] += monomial.coefficient;
                    }
                    degree = std::max(degree, monomial.degree);
                    if (monomial.term.size()) {
                        maxIndex = std::max(maxIndex, *std::max_element(monomial.term.begin(), monomial.term.end()));
                    }
                }
            }

            Polynomial(std::string fileName)
            {
                std::ifstream ifs;
                ifs.open(fileName);

                int termNum;
                ifs >> termNum;

                this->degree = 0;
                this->maxIndex = -1;

                for (int i = 0; i < termNum; i++) {
                    int deg;
                    ifs >> deg;
                    this->degree = std::max(this->degree, deg);

                    Term tmp(deg);
                    for (int j = 0; j < deg; j++) {
                        ifs >> tmp[j]; tmp[j]--;
                    }
                    std::sort(tmp.begin(), tmp.end());
                    if (tmp.size()) {
                        this->maxIndex = std::max(this->maxIndex, *std::max_element(tmp.begin(), tmp.end()));
                    }

                    double coef;
                    ifs >> coef;
                    if (monomials.find(tmp) == monomials.end()) {
                        monomials[tmp] = coef;
                    }
                    else {
                        monomials[tmp] += coef;
                    }
                }
            }

            void LoadFromFile(std::string fileName)
            {
                Polynomial tmp(fileName);
                *this = std::move(tmp); 
            }

            Polynomial operator + (const Polynomial &rhs)
            {
                Polynomial ret = *this;
                ret.degree = std::max(this->degree, rhs.degree);
                ret.maxIndex = std::max(this->maxIndex, rhs.maxIndex);
                for (auto monomial : rhs.monomials) {
                    if (ret.monomials.find(monomial.first) == monomials.end()) {
                        ret.monomials[monomial.first] = monomial.second;
                    }
                    else {
                        ret.monomials[monomial.first] += monomial.second;
                    }
                }
                return ret;
            }

            Polynomial &operator += (const Polynomial &rhs)
            {
                *this = *this + rhs;
                return *this;
            }

            const Polynomial operator -() const
            {
                Polynomial ret = *this;
                for (auto &monomial : ret.monomials) {
                    monomial.second = -monomial.second;
                }
                return ret;
            }

            Polynomial operator - (const Polynomial &rhs)
            {
                return *this + (-rhs);
            }

            Polynomial &operator -= (const Polynomial &rhs)
            {
                *this = *this - rhs;
                return *this;
            }

            Polynomial operator * (const Monomial &m)
            {
                Polynomial ret;
                for (auto monomial : monomials) {
                    Monomial newMonomial = Monomial(monomial.first, monomial.second) * m;
                    ret.monomials[newMonomial.term] = newMonomial.coefficient;
                }
                ret.degree = this->degree + m.degree;
                if (m.term.size()) {
                    ret.maxIndex = std::max(this->maxIndex, *std::max_element(m.term.begin(), m.term.end()));
                }
                return ret;
            }

            Polynomial operator * (const Polynomial &p)
            {
                Polynomial ret;
                for (auto monomial : p.monomials) {
                    Polynomial tmp = *this * Monomial(monomial.first, monomial.second);
                    ret = ret + tmp;
                }
                return ret;
            }

            Polynomial DifferentiateBy(int ind)
            {
                Polynomial ret = Monomial();
                for (auto &monomial : this->monomials) {
                    ret = ret + Monomial(monomial.first, monomial.second).DifferentiateBy(ind);
                }
                return ret;
            }

            std::string ToString(bool oneIndexed = false)
            {
                std::string res;
                int outputNum = 0;
                for (auto monomial : monomials) {
                    res += Monomial(monomial.first, monomial.second).ToString(/* oneIndexed = */oneIndexed, /* withSign = */outputNum != 0);
                    outputNum++;
                }
                return outputNum == 0 ? "0" : res;
            }

            double Evaluate(std::vector<double> &values)
            {
                double res = 0;
                for (auto &monomial : monomials) {
                    res += Monomial(monomial.first, monomial.second, /* sorted = */true).Evaluate(values);
                }
                return res;
            }

            int degree;
            int maxIndex;
            std::map<Term, double> monomials;
    };
}

#endif //__POLYNOMIAL_HPP__