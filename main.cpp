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
    BigInteger g("2410497055970432881345493397846112198995088771364307195189734031205605186951241875096796459061741781667380437076874705300974836586165283741668656930807264789");
    BigInteger ya("973768284341326272301553751114322685324340805902069613339667142187801529585352406975030927008752571917079716318221077758236884342662829402529734009607531649");
    BigInteger yb("4149822765985031146554298777122732051870868431387323913747785791685310508836836283702926446817000114868007555717546362425841865173929670156568682745060708314");
    BigInteger ord_g("4309874666");

    data d{p,g,ya,ord_g};

    clock_t t=clock();
    BigInteger xa=Pollard(d);
    t=clock()-t;
    xa.show();
    printf("Time elapse: %lf\n",1000.0*t/CLOCKS_PER_SEC);

    //    BigInteger temp=g;
//    temp=temp.mod(g);

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
