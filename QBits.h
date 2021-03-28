#ifndef __QBits_h__
#define __QBits_h__

#include<complex>
#include<vector>
#include<ostream>

#define NQBITS_MAX 22

typedef std::complex<double> cmplx;

struct QBits {// quantum computer simulator
    static cmplx *c;// pointer to amplitudes of basis states
    static long N;// number of currently used qubits
    static long D;// dimension (2^N)
    static long L;// number of allocated qubits
    long b;// begining index of this qubits
    long e;// ending index of this qubits + 1
    long n;// number of this qubits (e-b)
    long d;// number of states (2^n)
    long m;// mask ((d-1)<<b)
    static void resize(long);
    static void kill();
    QBits(long);
    ~QBits();
    inline long operator[](long i) const
    { return 1L<<(i<0?e+i:b+i); }
    // q[i] = bit flag of i-th qubit
    // negative index is allowed like python
    inline long operator()(long i, long j) const
    { return ((1L<<(i<0?e+i:b+i))-1)^((1L<<(j<0?e+j:b+j))-1); }
    // q(i,j) = bit flags of consecutive qubits from i to (j-1)th
    inline long operator()(long i) const { return i<<b; }
    // q(i) = bit flags of i-th basis state
    inline long begin() const { return b; }
    inline long end() const { return e; }
    inline long len() const { return n; }
    inline long dim() const { return d; }
    inline long mask() const { return m; }
};

void Hadamar(long, long=0);
void XNot(long, long=0);
void YRot(long, double, long=0);
void ZFlip(long);
void phase(long, double);
void swap(long, long, long=0);
long measure(long);
void set(const QBits&, long);

// q = target qubits, k = control mask
inline void Hadamar(const QBits& q, long k=0) { Hadamar(q.m, k); }
inline void XNot(const QBits& q, long k=0) { XNot(q.m, k); }
inline void YRot(const QBits& q, double a, long k=0) { YRot(q.m, a, k); }
inline void ZFlip(const QBits& q, long k=0) { ZFlip(k|q.m); }
inline void phase(const QBits& q, double a, long k=0) { phase(k|q.m, a); }
inline long measure(const QBits& q) { return measure(q.m)>>q.b; }
inline void clear(const QBits& q) { set(q,0); }

typedef std::vector<std::vector<cmplx> > mat_cmplx;
typedef std::vector<std::vector<double> > mat_double;

long amplitude(mat_cmplx&, const QBits&);
long probability(mat_double&, const QBits&);

std::ostream& operator<<(std::ostream&, const QBits&);

void QFT(const QBits&, bool=false);
inline void invQFT(const QBits& q) { QFT(q, true); }

long PhaseEst(void(const QBits&, long), const QBits&, long);

#endif // __QBits_h__