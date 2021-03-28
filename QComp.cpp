#include "QComp.h"
#include "tools.h"

// reference:
//   E. R. Johnston et al, "Programming Quantum Computers"
//   V. Vedral, A. Barenco and A. Ekert, quant-ph 951108

void operator+=(const QBits& x, const QBits& y)
// x+=y; ref. Johnston et al, fig5.8
{
    long i,j,m(MIN(x.len(), y.len()));
    for(i=0; i<m; i++)
        for(j=x.len()-1; j>=i; j--)
            XNot(x[j], x(i,j)|y[i]);
}

void operator-=(const QBits& x, const QBits& y)
// x-=y; ref. Johnston et al, p92
{
    long i,j,m(MIN(x.len(), y.len()));
    for(i=m-1; i>=0; i--)
        for(j=i; j<x.len(); j++)
            XNot(x[j], x(i,j)|y[i]);
}

void add(const QBits& x, long a, long k)
// x+=a; ref. Johnston et al, fig12.19
// k = flags of control qubits
{
    long i,j,l(bitlen(a));
    for(i=0; i<l; i++)
        if(a&(1L<<i))
            for(j=x.len()-1; j>=i; j--)
                XNot(x[j], x(i,j)|k);
}

void sub(const QBits& x, long a, long k)
// x-=a; ref. Johnston et al, fig12.19
// k = flags of control qubits
{
    long i,j;
    for(i=bitlen(a)-1; i>=0; i--)
        if(a&(1L<<i))
            for(j=i; j<x.len(); j++)
                XNot(x[j], x(i,j)|k);
}

void AddMod(const QBits& x, long a, long n, long k)
// x = x+a mod n (0<=x<n, 0<=a<n)
// k = flags of control qubits
// ref. Vedral, Barenco and Ekert, fig4
{
    long j(k|x[-1]);// sign bit of x
    QBits s(1);// assume s has been set to zero
    add(x,a,k);
    sub(x,n,k);
    XNot(s,j);// detect overflow
    add(x,n,k|s[0]);
    sub(x,a,k);
    XNot(s,j);
    XNot(s,k);// reset s to zero
    add(x,a,k);
}

void SubMod(const QBits& x, long a, long n, long k)
// x = x-a mod n (0<=x<n, 0<=a<n)
// k = flags of control qubits
// ref. Vedral, Barenco and Ekert, fig4
{
    long j(k|x[-1]);// sign bit of x
    QBits s(1);// assume s has been set to zero
    sub(x,a,k);
    XNot(s,k);
    XNot(s,j);// detect overflow
    add(x,a,k);
    sub(x,n,k|s[0]);
    XNot(s,j);// reset s to zero
    add(x,n,k);
    sub(x,a,k);
}

void MulAddMod(const QBits& y, const QBits& x, long a, long n, long k)
// y = y + x*a mod n (0<=x<n, 0<=a<n)
// k = flags of control qubits
// ref. Vedral, Barenco and Ekert, fig5
{
    long i;
    for(i=0; i<x.len(); i++) {
        AddMod(y,a,n,k|x[i]);
        a<<=1; a%=n;
    }
}

void MulSubMod(const QBits& y, const QBits& x, long a, long n, long k)
// y = y - x*a mod n (0<=x<n, 0<=a<n)
// k = flags of control qubits
// ref. Vedral, Barenco and Ekert, fig5
{
    long i;
    for(i=0; i<x.len(); i++) {
        SubMod(y,a,n,k|x[i]);
        a<<=1; a%=n;
    }
}

void swap(const QBits& x, const QBits& y, long k)
{
    long i,m(MIN(x.len(), y.len()));
    for(i=0; i<m; i++)
        swap(x[i], y[i], k);
}

void PowMulMod(const QBits& y, long a, const QBits& x, long n)
// y = y * a^x mod n
// ref. Vedral, Barenco and Ekert, fig6
{
    long i,j;
    QBits z(y.len()+1);// assume z has been set to zero
                       // extra 1 bit to detect overflow
    for(i=0; i<x.len(); i++) {
        j = x[i];
        MulAddMod(z,y,a,n,j);
        swap(y,z,j);// store result in y
        MulSubMod(z,y,InvMod(a,n),n,j);// reset z to zero
        a*=a; a%=n;
    }
}

void negate(const QBits& q) { XNot(q); q++; }