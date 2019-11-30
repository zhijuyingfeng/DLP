#include <ctime>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "bigint.h"

//Used in Pollard
struct tuple
{
    BigInteger x,a,b;
};

struct data
{
    BigInteger p,alpha,beta,n;
};

//Used in CRT
struct congruence
{
    BigInteger a,n;
};

//Used in SHANKS
struct map
{
    int i;
    BigInteger bg;
};

void f(tuple&t,const data&d);
BigInteger Pollard(const data&d);
BigInteger Pohlig_Hellman(const data&d,const BigInteger&q,const BigInteger&c);
BigInteger CRT(const congruence*c,const int&len);
BigInteger SHANKS(const data&d);
bool cmp(const map&a,const map&b);

int main()
{
    BigInteger p("5682549022748424631339131913370125786212509227588493537874673173634936008725904358935442101466555561124455782847468955028529037660533553941399408331331403379");
    BigInteger g("2410497055970432881345493397846112198995088771364307195189734031205605186951241875096796459061741781667380437076874705300974836586165283741668656930807264789");
    BigInteger ya("973768284341326272301553751114322685324340805902069613339667142187801529585352406975030927008752571917079716318221077758236884342662829402529734009607531649");
    BigInteger yb("4149822765985031146554298777122732051870868431387323913747785791685310508836836283702926446817000114868007555717546362425841865173929670156568682745060708314");
    BigInteger ord_g("4309874666");//2154937333*2

    data d{p,g,ya,ord_g};

    congruence c[2];
    c[0].n=BigInteger("2");
    c[1].n=BigInteger("2154937333");
    clock_t t=clock();
    c[0].a=Pohlig_Hellman(d,c[0].n,BigInteger::ONE);
    c[1].a=Pohlig_Hellman(d,c[1].n,BigInteger::ONE);
    BigInteger ans=CRT(c,2);
    ans=yb.modPow(ans,p);
    t=clock()-t;
    ans.show();
    printf("Time elapsed:\t%lf\n",1000.0*t/CLOCKS_PER_SEC);

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
    else
        return BigInteger::ZERO;
}

BigInteger Pohlig_Hellman(const data&d,const BigInteger&q,const BigInteger&c)
{
    BigInteger j(BigInteger::ZERO);
    BigInteger beta=d.beta;
    BigInteger temp_beta=d.n.divide(q);//n/(q^(j+1))
    BigInteger n_divide_q=d.n.divide(q);
    BigInteger ans(BigInteger::ZERO);
    BigInteger base(BigInteger::ONE);
    BigInteger gamma=d.alpha.modPow(n_divide_q,d.p);//alpha^(n/q)
    BigInteger alpha_inv=d.alpha.modInverse(d.p);
    BigInteger temp_q=alpha_inv;//alpha^(- q^j)
    while(!c.subtract(BigInteger::ONE).subtract(j).isNegative())
    {
        BigInteger sigma=beta.modPow(temp_beta,d.p);
        temp_beta=temp_beta.divide(q);
//        BigInteger gamma_pow_i(gamma);
//        for(BigInteger i=BigInteger::ONE;i.compareTo(q)<0;i=i.add(BigInteger::ONE))
//        {
//            if(gamma_pow_i.compareTo(sigma)==0)
//            {
//                ans=ans.add(i.multiply(base));
//                base=base.multiply(q);
//                beta=beta.multiply(d.alpha.modInverse(d.p).modPow(i.multiply(temp_q),d.p));
//                temp_q=temp_q.multiply(q);
//                break;
//            }
//            gamma_pow_i=gamma_pow_i.multiply(gamma).mod(d.p);
//        }
        data d_={d.p,gamma,sigma,q};
        BigInteger i=SHANKS(d_);
        ans=ans.add(base.multiply(i));
        base=base.multiply(q);
        beta=beta.multiply(temp_q.modPow(i,d.p));
        temp_q=temp_q.modPow(q,d.p);
        j=j.add(BigInteger::ONE);
    }
    return ans;
}

BigInteger CRT(const congruence*c,const int&len)
{
    BigInteger M=BigInteger::ONE;
    for(int i=0;i<len;i++)
        M=M.multiply(c[i].n);
    BigInteger *m=new BigInteger[len];
    BigInteger *t=new BigInteger[len];
    BigInteger ans(BigInteger::ZERO);
    for(int i=0;i<len;i++)
    {
        m[i]=M.divide(c[i].n);
        t[i]=m[i].modInverse(c[i].n);
        ans=ans.add(c[i].a.multiply(t[i]).multiply(m[i]));
    }
    delete [] m;
    delete [] t;
    return ans.mod(M);
}

BigInteger SHANKS(const data&d)
{
    int m=static_cast<int>(sqrt(d.n.longValue()))+1;
    map *L1=new map[m];
    map *L2=new map[m];
    char s[15]={0};
    sprintf(s,"%d",m);
    BigInteger m_(s);
    BigInteger alpha_m=d.alpha.modPow(m_,d.p);
    BigInteger alpha_inv=d.alpha.modInverse(d.p);
    BigInteger ans("-1");
    L1[0].bg=BigInteger::ONE;
    L1[0].i=0;
    L2[0].i=0;
    L2[0].bg=d.beta;
    for(int i=1;i<m;i++)
    {
        L1[i].bg=L1[i-1].bg.multiply(alpha_m).mod(d.p);
        L1[i].i=i;
        L2[i].i=i;
        L2[i].bg=L2[i-1].bg.multiply(alpha_inv).mod(d.p);
    }

    std::sort(L1,L1+m,cmp);
    std::sort(L2,L2+m,cmp);

    int i=0,j=0;
    while(i<m&&j<m)
    {
        int res=L1[i].bg.compareTo(L2[j].bg);
        if(res<0)
            i++;
        else if(res>0)
            j++;
        else
        {
            char str[30]={0};
            sprintf(str,"%lld",static_cast<long long>(m)*static_cast<long long>(L1[i].i)+static_cast<long long>(L2[j].i));
            ans=BigInteger(str).mod(d.n);
            break;
        }
    }
    delete [] L1;
    delete [] L2;
    return ans;
}

bool cmp(const map&a,const map&b)
{
    int res=a.bg.compareTo(b.bg);
    return res<0?1:0;
}
