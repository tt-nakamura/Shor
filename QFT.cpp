#include "QBits.h"
#include "tools.h"

void QFT(const QBits& q, bool inv)
// Quantum Fourier Transform
// inv = inverse QFT or not
// referece: Nielsen and Chuang, fig5.1
{
    long i,j,m(q[0]),n(q[-1]);
    double phi, pi2(inv ? PI/2:-PI/2);
    for(i=n; i>=m; i>>=1) {
        Hadamar(i);
        phi = pi2;
        for(j=i>>1; j>=m; j>>=1) {
            phase(i|j, phi);
            phi /= 2;
        }
    }
    while(m<n) {// reversing qubits
        swap(m,n);
        m<<=1; n>>=1;
    }
}

long PhaseEst(void U(const QBits&, long), const QBits& q, long t)
// phase estimation
// input:
//   U = unitary operator whose eigenphase is to be estimated
//   q = superposition of eigen state of U
//   t = number of qubits for phase estimation
// return:
//   eigen phase of U (in 2pi radian)
//     approximated as t bit integer
// reference: Nielsen and Chuang, fig5.2, 5.3
{
    long i,j,k(1);
    QBits q1(t);
    Hadamar(q1);
    for(i=0; i<t; i++, k<<=1)
        for(j=0; j<k; j++) U(q,q1[i]);
    invQFT(q1);
    i = measure(q1);
    XNot(q1(i));// clear q1 before destruct
    return i;
}