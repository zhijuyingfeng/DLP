#include <cstring>
#include "mpn.h"

int32_t MPN::NEGATIVE_ONE=-1;
int32_t MPN::CHARS_PER_WORD=9;
int64_t MPN::NEGATIVE_ONE_64=0xFFFFFFFF;

int32_t MPN::add_1(int32_t *dest, int32_t *x, const int32_t &size, const int32_t &y)
{
    int64_t carry=y;
    for(int i=0;i<size;i++)
    {
        carry+=static_cast<int64_t>(x[i])&NEGATIVE_ONE_64;
        dest[i]=carry&NEGATIVE_ONE;
        carry>>=32;
    }
    return carry&NEGATIVE_ONE;
}


int32_t MPN::mul_1(int32_t *dest, int32_t *x, const int32_t &len, const int32_t &y)
{
    int64_t carry=0;
    int64_t ywords=static_cast<int64_t>(y)&NEGATIVE_ONE_64;
    for(int i=0;i<len;i++)
    {
        carry+=(NEGATIVE_ONE_64&static_cast<int64_t>(x[i]))*ywords;
        dest[i]=carry&NEGATIVE_ONE;
        uint64_t temp=static_cast<uint64_t>(carry);
        temp>>=32;
        carry=static_cast<int64_t>(temp);
    }
    return NEGATIVE_ONE&carry;
}

int32_t MPN::set_str(int32_t *dest, int8_t *str, const int32_t &str_len)
{
    int32_t pos=0;
    int32_t base=10;
    int32_t size=0;

    while(pos<str_len)
    {
        int32_t temp=0;
        int32_t bigbase=1;
        for(int32_t i=0;i<CHARS_PER_WORD&&pos<str_len;i++,pos++)
        {
            temp=temp*base+str[pos];
            bigbase*=base;
        }
        int cy_limb;
        if(!size)
        {
            cy_limb=temp;
        }
        else
        {
            cy_limb=mul_1(dest,dest,size,bigbase);
            cy_limb+=add_1(dest,dest,size,temp);
        }
        if(cy_limb)
            dest[size++]=cy_limb;
    }
    return size;
}

int32_t MPN::add_n(int32_t *dest, int32_t *x, int32_t *y, const int32_t &len)
{
    int64_t carry=0;
    for(int i=0;i<len;i++)
    {
        carry+=(static_cast<int64_t>(x[i])&NEGATIVE_ONE_64)+
                (static_cast<int64_t>(y[i])&NEGATIVE_ONE_64);
        dest[i]=carry&NEGATIVE_ONE;
        uint64_t temp=static_cast<uint64_t>(carry);
        temp>>=32;
        carry=static_cast<int64_t>(temp);
    }
    return carry&NEGATIVE_ONE;
}

void MPN::mul(int32_t *dest, int32_t *x, const int32_t &xlen, int32_t *y, const int32_t &ylen)
{
    dest[xlen]=mul_1(dest,x,xlen,y[0]);

    for(int32_t i=1;i<ylen;i++)
    {
        int64_t carry=0,yword=static_cast<int64_t>(y[i])&NEGATIVE_ONE_64;

        for(int32_t j=0;j<xlen;j++)
        {
            int64_t temp=static_cast<int64_t>(x[j])&NEGATIVE_ONE_64;
            temp*=yword;
            temp+=static_cast<int64_t>(dest[i+j])&NEGATIVE_ONE_64;

            carry+=temp;
            dest[i+j]=carry&NEGATIVE_ONE;

            uint64_t tmp=static_cast<uint64_t>(carry);
            tmp>>=32;
            carry=static_cast<int64_t>(tmp);
        }
        dest[i+xlen]=carry&NEGATIVE_ONE;
    }
}

int32_t MPN::count_leading_zeros(const int32_t &x)
{
    if(x<0)
        return 0;
    uint32_t n=static_cast<uint32_t>(x);
    int32_t l=0,r=31,m;
    uint32_t temp;
    while(l<r)
    {
        m=(l+r)>>1;
        temp=n>>m;
        if(temp)
            l=m+1;
        else
            r=m;
    }
    return 32-l;
}

void MPN::normalize(int32_t *x, int32_t &xlen, int32_t *y,const int32_t &ylen)
{
    if(y[ylen-1]>=0)
    {
        int32_t y_zeros=count_leading_zeros(y[ylen-1]);
        int64_t overflow=0,temp;
        for(int32_t i=0;i<ylen;i++)
        {
            temp=static_cast<int64_t>(y[i])&NEGATIVE_ONE_64;
            temp<<=y_zeros;
            temp|=overflow;
            y[i]=temp&NEGATIVE_ONE;
            overflow=temp>>32;
        }
        overflow=0;
        for(int32_t i=0;i<xlen;i++)
        {
            temp=static_cast<int64_t>(x[i])&NEGATIVE_ONE_64;
            temp<<=y_zeros;
            temp|=overflow;
            x[i]=temp&NEGATIVE_ONE;
            overflow=temp>>32;
        }
        if(overflow)
            x[xlen++]=overflow&NEGATIVE_ONE;
    }
}

int32_t MPN::sub_n(int32_t *dest, const int32_t *X, const int32_t *Y, const int32_t &size)
{
//    int cy = 0;
//    for (int i = 0;  i < size;  i++)
//    {
//        int y = Y[i];
//        int x = X[i];
//        y += cy;	/* add previous carry to subtrahend */
//        // Invert the high-order bit, because: (unsigned) X > (unsigned) Y
//        // iff: (int) (X^0x80000000) > (int) (Y^0x80000000).
//        cy = (y^0x80000000) < (cy^0x80000000) ? 1 : 0;
//        y = x - y;
//        cy += (y^0x80000000) > (x ^ 0x80000000) ? 1 : 0;
//        dest[i] = y;
//    }
//    return cy;
    int64_t cy=0;
    for(int32_t i=0;i<size;i++)
    {
        int32_t x=X[i]&NEGATIVE_ONE_64;
        int32_t y=Y[i]&NEGATIVE_ONE_64;
        int64_t temp=x-y-cy;

        if(temp<0)
        {
            cy=1;
            temp+=0x100000000LL;
        }
        else
            cy=0;
        dest[i]=temp&NEGATIVE_ONE;
    }
    return cy&NEGATIVE_ONE;
}

void MPN::divide(int32_t *zds, const int32_t &nx, const int32_t *y, const int32_t &ny)
{
    int32_t j=nx;
    do
    {
        int32_t qhat;
        if(zds[j]==y[ny-1])
            qhat=-1;
        else
        {
            int64_t w=(static_cast<int64_t>(zds[j])<<32)+
                    (static_cast<int64_t>(zds[j-1])&NEGATIVE_ONE_64);
            qhat=static_cast<int32_t>(udiv_qrnnd(w,y[ny-1]));
        }
        if(qhat)
        {
            int32_t borrow=submul_1(zds,j-ny,y,ny,qhat);
            int32_t save=zds[j];
            int64_t num=(static_cast<int64_t>(save)&NEGATIVE_ONE_64)
                    -(static_cast<int64_t>(borrow)&NEGATIVE_ONE_64);
            while(num)
            {
                qhat--;
                int64_t carry=0;
                for(int32_t i=0;i<ny;i++)
                {
                    carry+=(static_cast<int64_t>(zds[j-ny+i])&NEGATIVE_ONE_64)
                            +(static_cast<int64_t>(y[i])&NEGATIVE_ONE_64);
                    zds[j-ny+i]=carry&NEGATIVE_ONE;

                    carry=logic_shift_right(carry,32);
                }
                zds[j]+=carry;
                num=carry-1;
            }
        }
        zds[j]=qhat;
    }while(--j>=ny);
}

int64_t MPN::udiv_qrnnd(int64_t N, int32_t D)
{
    int64_t q, r;
    int64_t a1 = logic_shift_right(N,32);
    int64_t a0 = N & NEGATIVE_ONE_64;
    if (D >= 0)
    {
        if (a1 < ((D - a1 - (logic_shift_right(a0,31))) & NEGATIVE_ONE_64))
        {
            /* dividend, divisor, and quotient are nonnegative */
            q = N / D;
            r = N % D;
        }
        else
        {
            /* Compute c1*2^32 + c0 = a1*2^32 + a0 - 2^31*d */
            int64_t c = N - (static_cast<int64_t>(D) << 31);
            /* Divide (c1*2^32 + c0) by d */
            q = c / D;
            r = c % D;
            /* Add 2^31 to quotient */
            q += static_cast<int64_t>(1)<< 31;
        }
      }
    else
    {
        uint32_t temp=static_cast<uint32_t>(D);
        temp>>=1;
        int64_t b1 = temp;	/* d/2, between 2^30 and 2^31 - 1 */
        //int64_t c1 = (a1 >> 1); /* A/2 */
        //int c0 = (a1 << 31) + (a0 >> 1);
        int64_t c = logic_shift_right(N,1);
        if (a1 < b1 || (a1 >> 1) < b1)
        {
            if (a1 < b1)
            {
                q = c / b1;
                r = c % b1;
            }
            else /* c1 < b1, so 2^31 <= (A/2)/b1 < 2^32 */
            {
                c = ~(c - (b1 << 32));
                q = c / b1;  /* (A/2) / (d/2) */
                r = c % b1;
                q = (~q) & NEGATIVE_ONE_64;    /* (A/2)/b1 */
                r = (b1 - 1) - r; /* r < b1 => new r >= 0 */
            }
            r = 2 * r + (a0 & 1);
            if ((D & 1) != 0)
            {
                if (r >= q)
                {
                    r = r - q;
                }
                else if (q - r <= (static_cast<int64_t>(D) & NEGATIVE_ONE_64))
                {
                    r = r - q + D;
                    q -= 1;
                }
                else
                {
                    r = r - q + D + D;
                    q -= 2;
                }
            }
        }
        else				/* Implies c1 = b1 */
        {				/* Hence a1 = d - 1 = 2*b1 - 1 */
            if (a0 >= (static_cast<int64_t>(-D) & NEGATIVE_ONE_64))
            {
                q = -1;
                r = a0 + D;
            }
            else
            {
                q = -2;
                r = a0 + D + D;
            }
        }
    }

    return (r << 32) | (q & 0xFFFFFFFFl);
}

int64_t MPN::logic_shift_right(const int64_t &N, const int32_t &r)
{
    uint64_t temp=static_cast<uint64_t>(N);
    temp>>=(r&0x3F);
    return static_cast<int64_t>(temp);
}

int32_t MPN::divmod_1(int32_t *quotient, int *dividend, const int &len, const int &divisor)
{
    int32_t i=len-1;
    int64_t r=dividend[i];
    if((r&NEGATIVE_ONE_64)>=(static_cast<int64_t>(divisor)&NEGATIVE_ONE_64))
        r=0;
    else
    {
        quotient[i--]=0;
        r<<=32;
    }
    for(;i>=0;i--)
    {
        int32_t n0=dividend[i];
        r=(r&~NEGATIVE_ONE_64)|(n0&NEGATIVE_ONE_64);
        r=udiv_qrnnd(r,divisor);
        quotient[i]=r&NEGATIVE_ONE;
    }
    return (r>>32)&NEGATIVE_ONE;
}


int32_t MPN::cmp(const int32_t *x, const int32_t *y, const int32_t &size)
{
    int32_t temp=size;
    while(--temp>=0)
    {
//        int32_t x_word=x[temp],y_word=y[temp];
//        if(x_word!=y_word)
//            return (x_word^0x80000000)>(y_word^0x8000000)?1:-1;

        uint32_t x_word=static_cast<uint32_t>(x[temp]);
        uint32_t y_word=static_cast<uint32_t>(y[temp]);
        if(x_word!=y_word)
            return x_word>y_word?1:-1;
    }
    return 0;
}

int32_t MPN::cmp(const int32_t *x, const int32_t &xlen, const int32_t *y, const int32_t &ylen)
{
    return xlen>ylen?1:xlen<ylen?-1:cmp(x,y,xlen);
}

int32_t MPN::lshift(int32_t *dest, int32_t d_offset, int32_t *x, const int32_t &len, const int32_t &count)
{
    int32_t count_2=32-count;
    int32_t i=len-1;
    int32_t high_word=x[i];
    uint32_t temp=static_cast<uint32_t>(high_word);
    temp>>=count_2;
    int32_t retval=static_cast<int32_t>(temp);
    d_offset++;
    while(--i>=0)
    {
        int32_t low_word=x[i];
        temp=static_cast<uint32_t>(low_word);
        temp>>=count_2;
        dest[d_offset+i]=(high_word<<count)|(static_cast<int32_t>(temp));
        high_word=low_word;
    }
    dest[d_offset+i]=high_word<<count;
    return retval;
}

void MPN::rshift0(int32_t *dest, int32_t *x, const int32_t &x_start, const int32_t &len, const int32_t &count)
{
    if(count>0)
        rshift(dest,x,x_start,len,count);
    else
    {
        for(int32_t i=0;i<len;i++)
            dest[i]=x[i+x_start];
    }
}

int32_t MPN::rshift(int32_t *dest, int32_t *x, const int32_t &x_start, const int32_t &len, const int32_t &count)
{
    int32_t count_2=32-count;
    int32_t low_word=x[x_start];
    int32_t retval=low_word<<count_2;
    int32_t i=1;
    for(;i<len;i++)
    {
        int32_t high_word=x[x_start+i];
        uint32_t temp=static_cast<uint32_t>(low_word);
        temp>>=count;
        dest[i-1]=(static_cast<int32_t>(temp))|(high_word<<count_2);
        low_word=high_word;
    }
    uint32_t temp=static_cast<uint32_t>(low_word);
    temp>>=count;
    dest[i-1]=static_cast<int32_t>(temp);
    return retval;
}

int32_t MPN::submul_1(int32_t *dest, const int32_t &offset, const int32_t *x, const int32_t &len, const int32_t &y)
{
    int64_t yl=static_cast<int64_t>(y)&NEGATIVE_ONE_64;
    const int32_t MIN_INT32=static_cast<int32_t>(0x80000000);
    int32_t carry=0,j=0;
    do
    {
        int64_t prod=yl*(static_cast<int64_t>(x[j])&NEGATIVE_ONE_64);
        int32_t prod_low=prod&NEGATIVE_ONE;
        int32_t prod_high=(prod>>32)&NEGATIVE_ONE;
        prod_low+=carry;
        // Invert the high-order bit, because: (unsigned) X > (unsigned) Y
        // iff: (int) (X^0x80000000) > (int) (Y^0x80000000).
        carry=((prod_low^MIN_INT32)<(carry^MIN_INT32)?1:0)+prod_high;

        int32_t x_j=dest[offset+j];
        prod_low=x_j-prod_low;
        if((prod_low^MIN_INT32)>(x_j^MIN_INT32))
            carry++;
        dest[offset+j]=prod_low;
    }while(++j<len);
    return carry;
}
