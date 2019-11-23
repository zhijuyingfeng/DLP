#include <ctime>
#include <cstring>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include "bigint.h"

struct tuple
{
    BigInteger x,a,b;
};

struct data
{
    BigInteger p,alpha,beta,n;
};

void f(tuple&t,const data&d);
BigInteger Pollard(const data&d);

uint64_t generateLong()
{
    int32_t a=rand(),b=rand();
    return (static_cast<uint64_t>(a)<<32)+static_cast<uint64_t>(b);
}

int main()
{
    BigInteger p("5682549022748424631339131913370125786212509227588493537874673173634936008725904358935442101466555561124455782847468955028529037660533553941399408331331403379");
    BigInteger g("1610199694236867902103084586156260484306820050663503013866639623097930249543168806526814507254407948055750960964971310568352400992321611017840341431973890026");
    BigInteger ya("2709055067184273601088156901754587645701655059414068318925571703597032854074470907589738912643523545972759650231168290651911762224455879680869674825283169042");
    BigInteger yb("2565110974102206303651459489762238044059734881714761356706955726474403236853736097565934108870316005404497531409553851279009315797458432849040636409910765205");
    BigInteger ord_g("1206102098588167138789");

//    BigInteger p("2147483647");
//    BigInteger g("16303069");
//    BigInteger ya("634053622");
//    BigInteger ord_g("119304647");

    data d{p,g,ya,ord_g};

//    BigInteger exp("28736");
//    BigInteger x=g.modPow(exp,p);
//    x.show();

    clock_t t=clock();
    BigInteger xa=Pollard(d);
    t=clock()-t;
    xa.show();
    printf("Time elapse: %lf\n",1000.0*t/CLOCKS_PER_SEC);
    return 0;
}

void f(tuple&t,const data&d)
{
    int32_t r=t.x.mod3();
    switch (r)
    {
    case 1:
        t.x=t.x.multiply(d.beta).mod(d.p);
        t.b=t.b.add(BigInteger::ONE).mod(d.n);
        break;
    case 0:
        t.x=t.x.multiply(t.x).mod(d.p);
        t.a=t.a.add(t.a).mod(d.n);
        t.b=t.b.add(t.b).mod(d.n);
        break;
    case 2:
        t.x=t.x.multiply(d.alpha).mod(d.p);
        t.a=t.a.add(BigInteger::ONE).mod(d.n);
    }
}

BigInteger Pollard(const data&d)
{
    tuple t{BigInteger::ONE,BigInteger::ZERO,BigInteger::ZERO};
    f(t,d);
    tuple t_=t;
    f(t_,d);
    while(t.x.compareTo(t_.x))
    {
        f(t,d);
        f(t_,d);
        f(t_,d);
    }
    BigInteger b_diff=t_.b.subtract(t.b).mod(d.n);
    BigInteger gcd=d.n.gcd(b_diff);
    BigInteger a_diff=t.a.subtract(t_.a).mod(d.n);
    if(gcd.isOne())
    {
        BigInteger inv=b_diff.modInverse(d.n);
        return a_diff.multiply(inv).mod(d.n);
    }
    else if(a_diff.mod(gcd).isZero())
    {
        BigInteger n=d.n.divide(gcd);
        BigInteger a=a_diff.divide(gcd);
        BigInteger b=b_diff.divide(gcd);
        BigInteger b_inv=b.modInverse(n);
        return a.multiply(b_inv).mod(n);
    }
    else
        return BigInteger::ZERO;
}
