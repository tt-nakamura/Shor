#include<cmath>
#include "QBits.h"
#include "tools.h"

cmplx *QBits::c(0);
long QBits::L(0);
long QBits::N(0);
long QBits::D;

void QBits::resize(long l)
// l = number of qubits
{
    if(l > NQBITS_MAX) error("l > NQBITS_MAX");
    if(l < L) error("l < L");// don't decrease
    L = l;
    if(c==0) {
        c = new cmplx[1L<<L];
        c[0] = 1.;
    }
    else {
        cmplx *t = new cmplx[1L<<L];
        for(long i=0; i<D; i++) t[i] = c[i];
        delete[] c;
        c = t;
    }
}

void QBits::kill() {
    delete[] c;
    c=0;
    L=N=0;
}

QBits::QBits(long k):
// get new k qubits
b(N), e(N+k), n(k), d(1L<<n), m((d-1)<<b)
{
    if((N = e) > L) resize(N);
    D = 1L<<N;
}

QBits::~QBits()
// push q back to unused qubits
{
    if(e==N) D = 1L<<(N=b);
}

void XNot(long j, long k)
// j = flags of target qubits
// k = flags of control qubits
{
    long i;
    for(i=0; i<QBits::D; i++)
        if((i&k)==k && (i^j) > i)
            SWAP(QBits::c[i^j], QBits::c[i]);
}

void Hadamar(long l, long k)
// l = flags of target qubits
// k = flags of control qubits
{
    long i,j,m,n;
    cmplx t;
    for(m=1; m<=l; m<<=1) {
        if((m&l)==0) continue;
        n = m<<1;
        for(i=n; i<=QBits::D; i+=n) {
            for(j=i-m; j<i; j++) {
                if((j&k)!=k) continue;
                t = QBits::c[j-m];
                QBits::c[j-m] = (t + QBits::c[j])*RSQRT2;
                QBits::c[j]   = (t - QBits::c[j])*RSQRT2;
            }
        }
    }
}

void ZFlip(long k)
// k = flags of control qubits
{
    long i;
    for(i=0; i<QBits::D; i++)
        if((i&k)==k) QBits::c[i] = -QBits::c[i];
}

void phase(long k, double a)
// k = flags of control qubits
// a = angle of rotation (radian) about z-axis
{
    long i;
    cmplx z(cos(a),sin(a));
    for(i=0; i<QBits::D; i++)
        if((i&k)==k) QBits::c[i] *= z;
}

void YRot(long l, double a, long k)
// l = flags of target qubits
// a = angle of rotation (radian) about y-axis
// k = flags of control qubits
{
    long i,j,m,n;
    double c(cos(a/2)), s(sin(a/2));
    cmplx t;
    for(m=1; m<=l; m<<=1) {
        if((m&l)==0) continue;
        n = m<<1;
        for(i=n; i<=QBits::D; i+=n) {
            for(j=i-m; j<i; j++) {
                if((j&k)!=k) continue;
                t = QBits::c[j-m];
                QBits::c[j-m] = c*t - s*QBits::c[j];
                QBits::c[j]   = s*t + c*QBits::c[j];
            }
        }
    }    
}

void swap(long j, long k, long l)
// j,k = flags of qubits to swap
// l = flags of control qubits
// assume weight(j)==1, weight(k)==1
{
    long i;
    for(i=0; i<QBits::D; i++)
        if((i&l)==l && i&j && (i&k)==0)
            SWAP(QBits::c[i], QBits::c[i^j^k]);
}

long measure(long k)
// k = flags of qubits to measure
// return: measurement result (bit flags)
// reference: N.D.Mermin,
//   "Quantum Computer Science" section 1.9
{
    long i,j;
    double s(0),r,p[QBits::D];
    for(i=0; i<QBits::D; i++) p[i] = 0.;
    for(i=0; i<QBits::D; i++)
        p[i&k] += norm(QBits::c[i]);
    r = frand();// throw dice
    for(i=0; i<QBits::D; i++)
        if((s += p[i]) > r) break;
    if(i==QBits::D) error("measurement failed");
    s = 1./sqrt(p[i]);
    for(j=0; j<QBits::D; j++)// normalize
        if((j&k)==i) QBits::c[j] *= s;
        else QBits::c[j] = 0.;

    return i;
}

void set(const QBits& q, long k)
// set pure state of value k
{
    long j;
    if((j = measure(q)) == k) return;
    XNot(q(j^k));
}

long amplitude(mat_cmplx& c, const QBits& q)
// c[i][j] = amplitude of j-th state in q (j=0, .., q.dim())
//           with other qubits in i-th state
// return k = number of states of other qubits (i=0,..,k)
//            excluding such i that c[i][j]==0 for all j
{
    long i,j,k(0),t,m(q[0]-1);
    long n(QBits::D/q.dim());
    for(i=0; i<n; i++) {
        t = ((i&~m)<<q.len())|i&m;
        for(j=0; j<q.dim(); j++)
            if(QBits::c[t|q(j)] != 0.) break;
        if(j==q.dim()) continue;// exclude if c==0
        c.resize(k+1);
        c[k].resize(q.dim());
        for(j=0; j<q.dim(); j++)
            c[k][j] = QBits::c[t|q(j)];
        k++;
    }
    return k;
}

long probability(mat_double& p, const QBits& q)
// p = |amplitude|^2, return same k as returned by amplitude(c,q)
{
    long i,j,k;
    mat_cmplx c;
    k = amplitude(c,q);
    p.resize(k);
    for(i=0; i<k; i++) p[i].resize(q.dim());
    for(i=0; i<k; i++)
        for(j=0; j<q.dim(); j++)
            p[i][j] = norm(c[i][j]);
    return k;
}


std::ostream& operator<<(std::ostream& s, const QBits& q)
// print amplitudes of basis states in q
{
    long i,j,k;
    mat_cmplx c;
    k = amplitude(c,q);
    for(i=0; i<k; i++) {
        s << i << ':';// index of states of other qubits
        for(j=0; j<q.dim(); j++)
            s << c[i][j];
        s << '\n';
    }
    return s;
}