#ifndef MPN_H
#define MPN_H

#include <cstdint>

class MPN
{
    static int32_t NEGATIVE_ONE;
    static int64_t NEGATIVE_ONE_64;
    static int32_t CHARS_PER_WORD;
public:
    static int32_t set_str(int32_t*dest, int8_t*str, const int32_t& str_len);

    /** Add x[0:size-1] and y, and write the size least
    * significant words of the result to dest.
    * Return carry, either 0 or 1.
    * All values are unsigned.
    * This is basically the same as gmp's mpn_add_1. */
    static int32_t add_1(int32_t *dest, int32_t*x, const int32_t&size, const int32_t&y);

    /** Multiply x[0:len-1] by y, and write the len least
    * significant words of the product to dest[0:len-1].
    * Return the most significant word of the product.
    * All values are treated as if they were unsigned
    * (i.e. masked with 0xffffffffL).
    * OK if dest==x (not sure if this is guaranteed for mpn_mul_1).
    * This function is basically the same as gmp's mpn_mul_1.
    */
    static int32_t mul_1(int32_t*dest, int32_t*x, const int32_t&len, const int32_t&y);

    /** Add x[0:len-1] and y[0:len-1] and write the len least
      * significant words of the result to dest[0:len-1].
      * All words are treated as unsigned.
      * Return the carry, either 0 or 1
      * This function is basically the same as gmp's mpn_add_n.*/
    static int32_t add_n(int32_t*dest,int32_t*x,int32_t*y,const int32_t &len);


    /**
     * Multiply x[0:xlen-1] and y[0:ylen-1], and
     * write the result to dest[0:xlen+ylen-1].
     * The destination has to have space for xlen+ylen words,
     * even if the result might be one limb smaller.
     * This function requires that xlen >= ylen.
     * The destination must be distinct from either input operands.
     * All operands are unsigned.
     * This function is basically the same gmp's mpn_mul. */
    static void mul(int32_t*dest,int32_t*x,const int32_t &xlen,int32_t*y,const int32_t &ylen);

    static int32_t count_leading_zeros(const int32_t&x);

    /*
     * x need an extra space
     * Assumes:  xlen>ylen.
     */
    static void normalize(int32_t*x,int32_t& xlen,int32_t*y,const int32_t&ylen);

    /** Divide zds[0:nx] by y[0:ny-1].
     * The remainder ends up in zds[0:ny-1].
     * The quotient ends up in zds[ny:nx].
     * Assumes:  nx>ny.
     * All operands are unsigned.
     */
    static void divide(int32_t* zds,const int32_t& nx,const int32_t* y,const int32_t& ny);


    /* Subtract x[0:len-1]*y from dest[offset:offset+len-1].
     * All values are treated as if unsigned.
     * @return the most significant word of
     * the product, minus borrow-out from the subtraction.
     */
    static int32_t submul_1(int32_t*dest,const int32_t&offset,const int32_t*x,const int32_t&len,const int32_t&y);


    /** Subtract Y[0:size-1] from X[0:size-1], and write
     * the size least significant words of the result to dest[0:size-1].
     * Return borrow, either 0 or 1.
     * This is basically the same as gmp's mpn_sub_n function.
     */
    static int32_t sub_n(int32_t*dest,const int32_t*X,const int32_t*Y,const int32_t&size);

    /* Divide (unsigned long) N by (unsigned int) D.
     * Returns (remainder << 32)+(unsigned int)(quotient).
     * Assumes (unsigned int)(N>>32) < (unsigned int)D.
     * Code transcribed from gmp-2.0's mpn_udiv_w_sdiv function.
     */
    static int64_t udiv_qrnnd (int64_t N, int32_t D);

    static int64_t logic_shift_right(const int64_t&N,const int32_t &r);


    /** Divide divident[0:len-1] by (unsigned int)divisor.
     * Write result into quotient[0:len-1.
     * Return the one-word (unsigned) remainder.
     * OK for quotient==dividend.
     */
    static int32_t divmod_1(int32_t* quotient,int* dividend,const int&len, const int&divisor);


    /** Compare x[0:size-1] with y[0:size-1], treating them as unsigned integers.
     * @result -1, 0, or 1 depending on if x<y, x==y, or x>y.
     * This is basically the same as gmp's mpn_cmp function.
     */
    static int32_t cmp(const int32_t*x,const int32_t*y,const int32_t&size);

    /** Compare x[0:xlen-1] with y[0:ylen-1], treating them as unsigned integers.
     * @result -1, 0, or 1 depending on if x<y, x==y, or x>y.
     */
    static int32_t cmp(const int32_t*x,const int32_t&xlen,const int32_t*y,const int32_t&ylen);

    /* Shift x[0:len-1] left by count bits, and store the len least
     * significant words of the result in dest[d_offset:d_offset+len-1].
     * Return the bits shifted out from the most significant digit.
     * Assumes 0 < count < 32.
     * OK if dest==x.
     */
    static int32_t lshift(int32_t*dest,int32_t d_offset,int32_t* x,const int32_t&len,const int32_t&count);

    /* Shift x[x_start:x_start+len-1] count bits to the "right"
     * (i.e. divide by 2**count).
     * Store the len least significant words of the result at dest.
     * OK if dest==x.
     * Assumes: 0 <= count < 32
     * Same as rshift, but handles count==0 (and has no return value).
     */
    static void rshift0(int32_t* dest,int32_t* x,const int32_t&x_start,const int32_t&len,const int32_t&count);

    /* Shift x[x_start:x_start+len-1] count bits to the "right"
     * (i.e. divide by 2**count).
     * Store the len least significant words of the result at dest.
     * The bits shifted out to the right are returned.
     * OK if dest==x.
     * Assumes: 0 < count < 32
     */
    static int32_t rshift(int32_t* dest,int32_t* x,const int32_t&x_start,const int32_t&len,const int32_t& count);
};
#endif // MPN_H
