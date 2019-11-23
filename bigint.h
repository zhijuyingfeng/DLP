#ifndef BIGINT_H
#define BIGINT_H

#include <cstdint>

class BigInteger
{
private:
    int32_t ival;
    int32_t *words;

    BigInteger(const int32_t&val);
    BigInteger(const int64_t&val);

    /** Destructively set the value of this to a long. */
    void set(const int64_t&y);
    void set(const BigInteger&y);
    void set(const int32_t *words, const int32_t &len);
    BigInteger canonicalize();
    void setInvert();
    void setAdd(const int32_t&y);
    void setShift(const BigInteger&x,const int32_t&count);
    void setShiftLeft(const BigInteger&x,int32_t count);
    void setShiftRight(const BigInteger&x,int32_t count);

    /** Calculate Greatest Common Divisor for non-negative ints. */
    static int32_t gcd(const int32_t&a,const int32_t&b);

    const static int32_t FLOOR=1;
    const static int32_t CEILING=2;
    const static int32_t TRUNCATE=3;
    const static int32_t ROUND=4;
    const static int32_t NEGATIVE_ONE=-1;
    const static int64_t NEGATIVE_ONE_64=0xFFFFFFFFLL;
    const static int32_t CHARS_PER_WORD=9;
    const static int32_t MIN_INT32=static_cast<int32_t>(0x80000000);
    const static int64_t MIN_INT64=static_cast<int64_t>(0x8000000000000000);

    static BigInteger valueOf(int8_t *digits, const int32_t &byte_len, const bool &negative);
    static int32_t wordsNeeded(const int32_t *words, const int32_t &len);
    static BigInteger valueOf(const int64_t& val);

    /** Set dest[0:len-1] to the negation of src[0:len-1].
    * Return true if overflow (i.e. if src is -2**(32*len-1)).
    * Ok for src==dest. */
    static bool negate(int32_t *dest, int32_t *src, const int32_t &len);

    /** Destructively set this to the negative of x.
    * It is OK if x==this.*/
    void setNegative(const BigInteger& x);

    /** Add two BigIntegers, yielding their sum as another BigInteger. */
    static BigInteger add(BigInteger x,BigInteger y,const int32_t&k);

    /** Make a canonicalized BigInteger from an array of words.
    * The array may be reused (without copying). */
    static BigInteger make(int32_t *words, const int32_t &len);

    /** Add a BigInteger and an int, yielding a new BigInteger. */
    static BigInteger add(const BigInteger&x,const int32_t&y);

    /** Add two ints, yielding a BigInteger. */
    static BigInteger add(const int32_t&x,const int32_t& y);

    /** Set this to the sum of x and y.
     * OK if x==this. */
    void setAdd(const BigInteger &x,const int32_t& y);

    /** Copy the abolute value of this into an array of words.
    * Assumes words.length >= (this.words == null ? 1 : this.ival).
    * Result is zero-extended, but need not be a valid 2's complement number.
    */
    void getAbsolute(int32_t* words)const;

    /** Return the logical (bit-wise) "and" of a BigInteger and an int. */
    static BigInteger And(const BigInteger&x,const int32_t&y);

    static void swap(BigInteger&x,BigInteger&y);
    static int32_t compareTo(const BigInteger&x,const BigInteger&y);

    static BigInteger neg(const BigInteger& x);
    static BigInteger times(const BigInteger& x,int32_t y);
    static BigInteger times(const BigInteger&x,const BigInteger& y);
    static void divide(int64_t x,int64_t y,BigInteger* quotient,BigInteger *remainder,const int32_t&rounding_mode);
    static BigInteger shift(const BigInteger&val,const int32_t&count);
    static BigInteger abs(const BigInteger&x);

    /** Divide two integers, yielding quotient and remainder.
    * @param x the numerator in the division
    * @param y the denominator in the division
    * @param quotient is set to the quotient of the result (iff quotient!=null)
    * @param remainder is set to the remainder of the result
    *  (iff remainder!=null)
    * @param rounding_mode one of FLOOR, CEILING, TRUNCATE, or ROUND.
    */
    static void divide(const BigInteger&x,const BigInteger&y,BigInteger* quotient,BigInteger *remainder,const int32_t&rounding_mode);
public:
    const static BigInteger ONE;
    const static BigInteger ZERO;
    const static BigInteger TEN;

    BigInteger();
    ~BigInteger();
    BigInteger(const BigInteger&val);
    BigInteger(const char* val);
    BigInteger operator=(const BigInteger&val);
    BigInteger add(const BigInteger&val)const;
    BigInteger subtract(const BigInteger& val)const;
    BigInteger multiply(const BigInteger&val)const;
    BigInteger divide(const BigInteger&val)const;
    BigInteger mod(const BigInteger&val)const;
    bool isNegative()const;
    void show();
    BigInteger modInverse(const BigInteger &val)const;
    BigInteger modPow(const BigInteger &exponent,const BigInteger& m)const;
    bool isZero()const;
    bool isOne()const;
    int64_t longValue()const;
    BigInteger negate();
    BigInteger And(const BigInteger&y)const;
    BigInteger shiftRight(const int32_t&n)const;
    BigInteger shiftLeft(const int32_t&n)const;
    BigInteger abs()const;

    //Assume *this and b is greater than b
    BigInteger gcd(const BigInteger& b)const;

    int32_t mod3()const;

    /** Destructively negate this. */
    void setNegative();
    int32_t compareTo(const BigInteger&val)const;
};

#endif // BIGINT_H
