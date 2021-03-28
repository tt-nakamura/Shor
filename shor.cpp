#include "QComp.h"
#include "tools.h"

long shor(long& d, long n)
// factor n by shor method.
// can factor n=15; can barely factor n=21;
// when n>21, the number of qubits is so large
//   that the simulator just won't work.
// assume n>0 is composite and square-free.
// output: d = factor of n (1<d<n)
// return 0 if d was found successfully by shor
// return 1 if d was found but not by shor
// return -1 if d was not found
// reference: N. D. Mermin
//   "Quantum Computer Science" section 3.7
{
    long l(bitlen(n));
    long m(n==15 ? l : l<<1);// Mermin, p82
    long M(1L<<m),a,s,t;
    QBits x(m),y(l);
    for(a=2; a<n; a++) {
        if((d = gcd(a,n)) > 1) return 1;
        clear(x);
        set(y,1);// y=1
        Hadamar(x);// x=0,1,..,M-1 superposed
        PowMulMod(y,a,x,n);// y*=b^x mod n
        QFT(x);
        d = measure(x);
        if(d==0) continue;
        ContFrac(s,t,d,M,n);
        if(PowerMod(a,t,n)!=1) continue;
        while((t&1)==0) {
            t >>= 1;// t is an even order of a mod n
            d = gcd(n, PowerMod(a,t,n)+1);
            if(d>1 && d<n) return 0;
        }
    }
    return -1;
}