#include "nxbase_core.h"
#include "skblas.h"
namespace blas
{
f77integer (* ssyr2)(f77char *uplo, f77integer *n, f77real *alpha, f77real *x, f77integer *incx, f77real *y, f77integer *incy, f77real *a, f77integer *lda) = &SSYR2;
f77integer (* drotg)(f77doublereal *da, f77doublereal *db, f77doublereal *c__, f77doublereal *s) = &DROTG;
f77integer (* cher)(f77char *uplo, f77integer *n, f77real *alpha, f77complex *x, f77integer *incx, f77complex *a, f77integer *lda) = &CHER;
f77integer (* zhemv)(f77char *uplo, f77integer *n, f77doublecomplex *alpha, f77doublecomplex *a, f77integer *lda, f77doublecomplex *x, f77integer *incx, f77doublecomplex *beta, f77doublecomplex *y, f77integer *incy) = &ZHEMV;
f77real (* snrm2)(f77integer *n, f77real *x, f77integer *incx) = &SNRM2;
f77integer (* zhemm)(f77char *side, f77char *uplo, f77integer *m, f77integer *n, f77doublecomplex *alpha, f77doublecomplex *a, f77integer *lda, f77doublecomplex *b, f77integer *ldb, f77doublecomplex *beta, f77doublecomplex *c__, f77integer *ldc) = &ZHEMM;
f77integer (* ztrmm)(f77char *side, f77char *uplo, f77char *transa, f77char *diag, f77integer *m, f77integer *n, f77doublecomplex *alpha, f77doublecomplex *a, f77integer *lda, f77doublecomplex *b, f77integer *ldb) = &ZTRMM;
f77doublereal (* ddot)(f77integer *n, f77doublereal *dx, f77integer *incx, f77doublereal *dy, f77integer *incy) = &DDOT;
f77integer (* ztrmv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77doublecomplex *a, f77integer *lda, f77doublecomplex *x, f77integer *incx) = &ZTRMV;
f77integer (* saxpy)(f77integer *n, f77real *sa, f77real *sx, f77integer *incx, f77real *sy, f77integer *incy) = &SAXPY;
f77integer (* cher2k)(f77char *uplo, f77char *trans, f77integer *n, f77integer *k, f77complex *alpha, f77complex *a, f77integer *lda, f77complex *b, f77integer *ldb, f77real *beta, f77complex *c__, f77integer *ldc) = &CHER2K;
f77real (* sasum)(f77integer *n, f77real *sx, f77integer *incx) = &SASUM;
f77integer (* isamax)(f77integer *n, f77real *sx, f77integer *incx) = &ISAMAX;
f77integer (* csyr2k)(f77char *uplo, f77char *trans, f77integer *n, f77integer *k, f77complex *alpha, f77complex *a, f77integer *lda, f77complex *b, f77integer *ldb, f77complex *beta, f77complex *c__, f77integer *ldc) = &CSYR2K;
f77integer (* caxpy)(f77integer *n, f77complex *ca, f77complex *cx, f77integer *incx, f77complex *cy, f77integer *incy) = &CAXPY;
f77integer (* cher2)(f77char *uplo, f77integer *n, f77complex *alpha, f77complex *x, f77integer *incx, f77complex *y, f77integer *incy, f77complex *a, f77integer *lda) = &CHER2;
f77integer (* dger)(f77integer *m, f77integer *n, f77doublereal *alpha, f77doublereal *x, f77integer *incx, f77doublereal *y, f77integer *incy, f77doublereal *a, f77integer *lda) = &DGER;
f77integer (* sgemv)(f77char *trans, f77integer *m, f77integer *n, f77real *alpha, f77real *a, f77integer *lda, f77real *x, f77integer *incx, f77real *beta, f77real *y, f77integer *incy) = &SGEMV;
f77integer (* scopy)(f77integer *n, f77real *sx, f77integer *incx, f77real *sy, f77integer *incy) = &SCOPY;
f77integer (* csyrk)(f77char *uplo, f77char *trans, f77integer *n, f77integer *k, f77complex *alpha, f77complex *a, f77integer *lda, f77complex *beta, f77complex *c__, f77integer *ldc) = &CSYRK;
f77integer (* chemv)(f77char *uplo, f77integer *n, f77complex *alpha, f77complex *a, f77integer *lda, f77complex *x, f77integer *incx, f77complex *beta, f77complex *y, f77integer *incy) = &CHEMV;
f77integer (* zrotg)(f77doublecomplex *ca, f77doublecomplex *cb, f77doublereal *c__, f77doublecomplex *s) = &ZROTG;
f77integer (* csymm)(f77char *side, f77char *uplo, f77integer *m, f77integer *n, f77complex *alpha, f77complex *a, f77integer *lda, f77complex *b, f77integer *ldb, f77complex *beta, f77complex *c__, f77integer *ldc) = &CSYMM;
f77integer (* chemm)(f77char *side, f77char *uplo, f77integer *m, f77integer *n, f77complex *alpha, f77complex *a, f77integer *lda, f77complex *b, f77integer *ldb, f77complex *beta, f77complex *c__, f77integer *ldc) = &CHEMM;
f77integer (* sgemm)(f77char *transa, f77char *transb, f77integer *m, f77integer *n, f77integer *k, f77real *alpha, f77real *a, f77integer *lda, f77real *b, f77integer *ldb, f77real *beta, f77real *c__, f77integer *ldc) = &SGEMM;
f77integer (* zhbmv)(f77char *uplo, f77integer *n, f77integer *k, f77doublecomplex *alpha, f77doublecomplex *a, f77integer *lda, f77doublecomplex *x, f77integer *incx, f77doublecomplex *beta, f77doublecomplex *y, f77integer *incy) = &ZHBMV;
f77integer (* dsbmv)(f77char *uplo, f77integer *n, f77integer *k, f77doublereal *alpha, f77doublereal *a, f77integer *lda, f77doublereal *x, f77integer *incx, f77doublereal *beta, f77doublereal *y, f77integer *incy) = &DSBMV;
f77doublereal (* dnrm2)(f77integer *n, f77doublereal *x, f77integer *incx) = &DNRM2;
f77real (* scnrm2)(f77integer *n, f77complex *x, f77integer *incx) = &SCNRM2;
f77integer (* ztbmv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77integer *k, f77doublecomplex *a, f77integer *lda, f77doublecomplex *x, f77integer *incx) = &ZTBMV;
f77integer (* dspmv)(f77char *uplo, f77integer *n, f77doublereal *alpha, f77doublereal *ap, f77doublereal *x, f77integer *incx, f77doublereal *beta, f77doublereal *y, f77integer *incy) = &DSPMV;
f77integer (* zaxpy)(f77integer *n, f77doublecomplex *za, f77doublecomplex *zx, f77integer *incx, f77doublecomplex *zy, f77integer *incy) = &ZAXPY;
f77integer (* zher2)(f77char *uplo, f77integer *n, f77doublecomplex *alpha, f77doublecomplex *x, f77integer *incx, f77doublecomplex *y, f77integer *incy, f77doublecomplex *a, f77integer *lda) = &ZHER2;
f77integer (* dtbsv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77integer *k, f77doublereal *a, f77integer *lda, f77doublereal *x, f77integer *incx) = &DTBSV;
f77integer (* ssyrk)(f77char *uplo, f77char *trans, f77integer *n, f77integer *k, f77real *alpha, f77real *a, f77integer *lda, f77real *beta, f77real *c__, f77integer *ldc) = &SSYRK;
f77integer (* dgbmv)(f77char *trans, f77integer *m, f77integer *n, f77integer *kl, f77integer *ku, f77doublereal *alpha, f77doublereal *a, f77integer *lda, f77doublereal *x, f77integer *incx, f77doublereal *beta, f77doublereal *y, f77integer *incy) = &DGBMV;
f77integer (* idamax)(f77integer *n, f77doublereal *dx, f77integer *incx) = &IDAMAX;
f77integer (* dtrmm)(f77char *side, f77char *uplo, f77char *transa, f77char *diag, f77integer *m, f77integer *n, f77doublereal *alpha, f77doublereal *a, f77integer *lda, f77doublereal *b, f77integer *ldb) = &DTRMM;
f77integer (* zgeru)(f77integer *m, f77integer *n, f77doublecomplex *alpha, f77doublecomplex *x, f77integer *incx, f77doublecomplex *y, f77integer *incy, f77doublecomplex *a, f77integer *lda) = &ZGERU;
f77integer (* cgbmv)(f77char *trans, f77integer *m, f77integer *n, f77integer *kl, f77integer *ku, f77complex *alpha, f77complex *a, f77integer *lda, f77complex *x, f77integer *incx, f77complex *beta, f77complex *y, f77integer *incy) = &CGBMV;
f77integer (* zgemm)(f77char *transa, f77char *transb, f77integer *m, f77integer *n, f77integer *k, f77doublecomplex *alpha, f77doublecomplex *a, f77integer *lda, f77doublecomplex *b, f77integer *ldb, f77doublecomplex *beta, f77doublecomplex *c__, f77integer *ldc) = &ZGEMM;
f77integer (* chpmv)(f77char *uplo, f77integer *n, f77complex *alpha, f77complex *ap, f77complex *x, f77integer *incx, f77complex *beta, f77complex *y, f77integer *incy) = &CHPMV;
f77integer (* stpsv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77real *ap, f77real *x, f77integer *incx) = &STPSV;
f77integer (* dsyrk)(f77char *uplo, f77char *trans, f77integer *n, f77integer *k, f77doublereal *alpha, f77doublereal *a, f77integer *lda, f77doublereal *beta, f77doublereal *c__, f77integer *ldc) = &DSYRK;
f77integer (* stbmv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77integer *k, f77real *a, f77integer *lda, f77real *x, f77integer *incx) = &STBMV;
f77integer (* zgemv)(f77char *trans, f77integer *m, f77integer *n, f77doublecomplex *alpha, f77doublecomplex *a, f77integer *lda, f77doublecomplex *x, f77integer *incx, f77doublecomplex *beta, f77doublecomplex *y, f77integer *incy) = &ZGEMV;
f77integer (* dtrmv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77doublereal *a, f77integer *lda, f77doublereal *x, f77integer *incx) = &DTRMV;
f77integer (* zgerc)(f77integer *m, f77integer *n, f77doublecomplex *alpha, f77doublecomplex *x, f77integer *incx, f77doublecomplex *y, f77integer *incy, f77doublecomplex *a, f77integer *lda) = &ZGERC;
f77integer (* sspr)(f77char *uplo, f77integer *n, f77real *alpha, f77real *x, f77integer *incx, f77real *ap) = &SSPR;
f77integer (* sscal)(f77integer *n, f77real *sa, f77real *sx, f77integer *incx) = &SSCAL;
f77integer (* ztpsv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77doublecomplex *ap, f77doublecomplex *x, f77integer *incx) = &ZTPSV;
f77integer (* chpr2)(f77char *uplo, f77integer *n, f77complex *alpha, f77complex *x, f77integer *incx, f77complex *y, f77integer *incy, f77complex *ap) = &CHPR2;
f77integer (* zsymm)(f77char *side, f77char *uplo, f77integer *m, f77integer *n, f77doublecomplex *alpha, f77doublecomplex *a, f77integer *lda, f77doublecomplex *b, f77integer *ldb, f77doublecomplex *beta, f77doublecomplex *c__, f77integer *ldc) = &ZSYMM;
f77integer (* zsyr2k)(f77char *uplo, f77char *trans, f77integer *n, f77integer *k, f77doublecomplex *alpha, f77doublecomplex *a, f77integer *lda, f77doublecomplex *b, f77integer *ldb, f77doublecomplex *beta, f77doublecomplex *c__, f77integer *ldc) = &ZSYR2K;
f77integer (* zhpr2)(f77char *uplo, f77integer *n, f77doublecomplex *alpha, f77doublecomplex *x, f77integer *incx, f77doublecomplex *y, f77integer *incy, f77doublecomplex *ap) = &ZHPR2;
f77integer (* ctrmm)(f77char *side, f77char *uplo, f77char *transa, f77char *diag, f77integer *m, f77integer *n, f77complex *alpha, f77complex *a, f77integer *lda, f77complex *b, f77integer *ldb) = &CTRMM;
f77integer (* ssyr)(f77char *uplo, f77integer *n, f77real *alpha, f77real *x, f77integer *incx, f77real *a, f77integer *lda) = &SSYR;
f77integer (* ctrmv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77complex *a, f77integer *lda, f77complex *x, f77integer *incx) = &CTRMV;
f77integer (* zswap)(f77integer *n, f77doublecomplex *zx, f77integer *incx, f77doublecomplex *zy, f77integer *incy) = &ZSWAP;
f77integer (* icamax)(f77integer *n, f77complex *cx, f77integer *incx) = &ICAMAX;
f77integer (* ssymm)(f77char *side, f77char *uplo, f77integer *m, f77integer *n, f77real *alpha, f77real *a, f77integer *lda, f77real *b, f77integer *ldb, f77real *beta, f77real *c__, f77integer *ldc) = &SSYMM;
f77integer (* ssymv)(f77char *uplo, f77integer *n, f77real *alpha, f77real *a, f77integer *lda, f77real *x, f77integer *incx, f77real *beta, f77real *y, f77integer *incy) = &SSYMV;
f77integer (* crotg)(f77complex *ca, f77complex *cb, f77real *c__, f77complex *s) = &CROTG;
f77integer (* dgemm)(f77char *transa, f77char *transb, f77integer *m, f77integer *n, f77integer *k, f77doublereal *alpha, f77doublereal *a, f77integer *lda, f77doublereal *b, f77integer *ldb, f77doublereal *beta, f77doublereal *c__, f77integer *ldc) = &DGEMM;
f77integer (* dsyr2)(f77char *uplo, f77integer *n, f77doublereal *alpha, f77doublereal *x, f77integer *incx, f77doublereal *y, f77integer *incy, f77doublereal *a, f77integer *lda) = &DSYR2;
f77integer (* dgemv)(f77char *trans, f77integer *m, f77integer *n, f77doublereal *alpha, f77doublereal *a, f77integer *lda, f77doublereal *x, f77integer *incx, f77doublereal *beta, f77doublereal *y, f77integer *incy) = &DGEMV;
f77integer (* zhpmv)(f77char *uplo, f77integer *n, f77doublecomplex *alpha, f77doublecomplex *ap, f77doublecomplex *x, f77integer *incx, f77doublecomplex *beta, f77doublecomplex *y, f77integer *incy) = &ZHPMV;
f77integer (* zher)(f77char *uplo, f77integer *n, f77doublereal *alpha, f77doublecomplex *x, f77integer *incx, f77doublecomplex *a, f77integer *lda) = &ZHER;
f77integer (* dsyr2k)(f77char *uplo, f77char *trans, f77integer *n, f77integer *k, f77doublereal *alpha, f77doublereal *a, f77integer *lda, f77doublereal *b, f77integer *ldb, f77doublereal *beta, f77doublereal *c__, f77integer *ldc) = &DSYR2K;
f77doublereal (* dznrm2)(f77integer *n, f77doublecomplex *x, f77integer *incx) = &DZNRM2;
f77integer (* sspr2)(f77char *uplo, f77integer *n, f77real *alpha, f77real *x, f77integer *incx, f77real *y, f77integer *incy, f77real *ap) = &SSPR2;
f77integer (* strmv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77real *a, f77integer *lda, f77real *x, f77integer *incx) = &STRMV;
f77integer (* ctrsv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77complex *a, f77integer *lda, f77complex *x, f77integer *incx) = &CTRSV;
f77integer (* dsymm)(f77char *side, f77char *uplo, f77integer *m, f77integer *n, f77doublereal *alpha, f77doublereal *a, f77integer *lda, f77doublereal *b, f77integer *ldb, f77doublereal *beta, f77doublereal *c__, f77integer *ldc) = &DSYMM;
f77integer (* zdscal)(f77integer *n, f77doublereal *da, f77doublecomplex *zx, f77integer *incx) = &ZDSCAL;
f77integer (* dsymv)(f77char *uplo, f77integer *n, f77doublereal *alpha, f77doublereal *a, f77integer *lda, f77doublereal *x, f77integer *incx, f77doublereal *beta, f77doublereal *y, f77integer *incy) = &DSYMV;
f77integer (* ztpmv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77doublecomplex *ap, f77doublecomplex *x, f77integer *incx) = &ZTPMV;
f77integer (* dsyr)(f77char *uplo, f77integer *n, f77doublereal *alpha, f77doublereal *x, f77integer *incx, f77doublereal *a, f77integer *lda) = &DSYR;
f77integer (* daxpy)(f77integer *n, f77doublereal *da, f77doublereal *dx, f77integer *incx, f77doublereal *dy, f77integer *incy) = &DAXPY;
f77integer (* ctrsm)(f77char *side, f77char *uplo, f77char *transa, f77char *diag, f77integer *m, f77integer *n, f77complex *alpha, f77complex *a, f77integer *lda, f77complex *b, f77integer *ldb) = &CTRSM;
f77integer (* strmm)(f77char *side, f77char *uplo, f77char *transa, f77char *diag, f77integer *m, f77integer *n, f77real *alpha, f77real *a, f77integer *lda, f77real *b, f77integer *ldb) = &STRMM;
f77integer (* zsyrk)(f77char *uplo, f77char *trans, f77integer *n, f77integer *k, f77doublecomplex *alpha, f77doublecomplex *a, f77integer *lda, f77doublecomplex *beta, f77doublecomplex *c__, f77integer *ldc) = &ZSYRK;
f77integer (* dspr2)(f77char *uplo, f77integer *n, f77doublereal *alpha, f77doublereal *x, f77integer *incx, f77doublereal *y, f77integer *incy, f77doublereal *ap) = &DSPR2;
f77integer (* drot)(f77integer *n, f77doublereal *dx, f77integer *incx, f77doublereal *dy, f77integer *incy, f77doublereal *c__, f77doublereal *s) = &DROT;
f77integer (* zherk)(f77char *uplo, f77char *trans, f77integer *n, f77integer *k, f77doublereal *alpha, f77doublecomplex *a, f77integer *lda, f77doublereal *beta, f77doublecomplex *c__, f77integer *ldc) = &ZHERK;
f77integer (* ssbmv)(f77char *uplo, f77integer *n, f77integer *k, f77real *alpha, f77real *a, f77integer *lda, f77real *x, f77integer *incx, f77real *beta, f77real *y, f77integer *incy) = &SSBMV;
f77integer (* chbmv)(f77char *uplo, f77integer *n, f77integer *k, f77complex *alpha, f77complex *a, f77integer *lda, f77complex *x, f77integer *incx, f77complex *beta, f77complex *y, f77integer *incy) = &CHBMV;
f77integer (* srot)(f77integer *n, f77real *sx, f77integer *incx, f77real *sy, f77integer *incy, f77real *c__, f77real *s) = &SROT;
f77integer (* sgbmv)(f77char *trans, f77integer *m, f77integer *n, f77integer *kl, f77integer *ku, f77real *alpha, f77real *a, f77integer *lda, f77real *x, f77integer *incx, f77real *beta, f77real *y, f77integer *incy) = &SGBMV;
f77integer (* izamax)(f77integer *n, f77doublecomplex *zx, f77integer *incx) = &IZAMAX;
f77integer (* ctpsv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77complex *ap, f77complex *x, f77integer *incx) = &CTPSV;
f77integer (* ctbmv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77integer *k, f77complex *a, f77integer *lda, f77complex *x, f77integer *incx) = &CTBMV;
f77integer (* cgemm)(f77char *transa, f77char *transb, f77integer *m, f77integer *n, f77integer *k, f77complex *alpha, f77complex *a, f77integer *lda, f77complex *b, f77integer *ldb, f77complex *beta, f77complex *c__, f77integer *ldc) = &CGEMM;
f77integer (* stpmv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77real *ap, f77real *x, f77integer *incx) = &STPMV;
f77integer (* cgeru)(f77integer *m, f77integer *n, f77complex *alpha, f77complex *x, f77integer *incx, f77complex *y, f77integer *incy, f77complex *a, f77integer *lda) = &CGERU;
f77integer (* zhpr)(f77char *uplo, f77integer *n, f77doublereal *alpha, f77doublecomplex *x, f77integer *incx, f77doublecomplex *ap) = &ZHPR;
f77integer (* cgemv)(f77char *trans, f77integer *m, f77integer *n, f77complex *alpha, f77complex *a, f77integer *lda, f77complex *x, f77integer *incx, f77complex *beta, f77complex *y, f77integer *incy) = &CGEMV;
f77doublereal (* dzasum)(f77integer *n, f77doublecomplex *zx, f77integer *incx) = &DZASUM;
f77integer (* cgerc)(f77integer *m, f77integer *n, f77complex *alpha, f77complex *x, f77integer *incx, f77complex *y, f77integer *incy, f77complex *a, f77integer *lda) = &CGERC;
f77integer (* zher2k)(f77char *uplo, f77char *trans, f77integer *n, f77integer *k, f77doublecomplex *alpha, f77doublecomplex *a, f77integer *lda, f77doublecomplex *b, f77integer *ldb, f77doublereal *beta, f77doublecomplex *c__, f77integer *ldc) = &ZHER2K;
f77integer (* ccopy)(f77integer *n, f77complex *cx, f77integer *incx, f77complex *cy, f77integer *incy) = &CCOPY;
f77integer (* srotg)(f77real *sa, f77real *sb, f77real *c__, f77real *s) = &SROTG;
f77integer (* cscal)(f77integer *n, f77complex *ca, f77complex *cx, f77integer *incx) = &CSCAL;
f77integer (* csscal)(f77integer *n, f77real *sa, f77complex *cx, f77integer *incx) = &CSSCAL;
f77integer (* zscal)(f77integer *n, f77doublecomplex *za, f77doublecomplex *zx, f77integer *incx) = &ZSCAL;
f77integer (* sspmv)(f77char *uplo, f77integer *n, f77real *alpha, f77real *ap, f77real *x, f77integer *incx, f77real *beta, f77real *y, f77integer *incy) = &SSPMV;
f77integer (* zgbmv)(f77char *trans, f77integer *m, f77integer *n, f77integer *kl, f77integer *ku, f77doublecomplex *alpha, f77doublecomplex *a, f77integer *lda, f77doublecomplex *x, f77integer *incx, f77doublecomplex *beta, f77doublecomplex *y, f77integer *incy) = &ZGBMV;
f77integer (* strsm)(f77char *side, f77char *uplo, f77char *transa, f77char *diag, f77integer *m, f77integer *n, f77real *alpha, f77real *a, f77integer *lda, f77real *b, f77integer *ldb) = &STRSM;
f77integer (* ssyr2k)(f77char *uplo, f77char *trans, f77integer *n, f77integer *k, f77real *alpha, f77real *a, f77integer *lda, f77real *b, f77integer *ldb, f77real *beta, f77real *c__, f77integer *ldc) = &SSYR2K;
f77integer (* strsv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77real *a, f77integer *lda, f77real *x, f77integer *incx) = &STRSV;
f77integer (* dscal)(f77integer *n, f77doublereal *da, f77doublereal *dx, f77integer *incx) = &DSCAL;
f77integer (* ctpmv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77complex *ap, f77complex *x, f77integer *incx) = &CTPMV;
f77integer (* dtpmv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77doublereal *ap, f77doublereal *x, f77integer *incx) = &DTPMV;
f77real (* scasum)(f77integer *n, f77complex *cx, f77integer *incx) = &SCASUM;
f77integer (* ztbsv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77integer *k, f77doublecomplex *a, f77integer *lda, f77doublecomplex *x, f77integer *incx) = &ZTBSV;
f77real (* sdot)(f77integer *n, f77real *sx, f77integer *incx, f77real *sy, f77integer *incy) = &SDOT;
f77integer (* dtpsv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77doublereal *ap, f77doublereal *x, f77integer *incx) = &DTPSV;
f77integer (* cswap)(f77integer *n, f77complex *cx, f77integer *incx, f77complex *cy, f77integer *incy) = &CSWAP;
f77integer (* dtrsv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77doublereal *a, f77integer *lda, f77doublereal *x, f77integer *incx) = &DTRSV;
f77integer (* ztrsv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77doublecomplex *a, f77integer *lda, f77doublecomplex *x, f77integer *incx) = &ZTRSV;
f77integer (* stbsv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77integer *k, f77real *a, f77integer *lda, f77real *x, f77integer *incx) = &STBSV;
f77integer (* dswap)(f77integer *n, f77doublereal *dx, f77integer *incx, f77doublereal *dy, f77integer *incy) = &DSWAP;
f77integer (* cherk)(f77char *uplo, f77char *trans, f77integer *n, f77integer *k, f77real *alpha, f77complex *a, f77integer *lda, f77real *beta, f77complex *c__, f77integer *ldc) = &CHERK;
f77integer (* ctbsv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77integer *k, f77complex *a, f77integer *lda, f77complex *x, f77integer *incx) = &CTBSV;
f77integer (* ztrsm)(f77char *side, f77char *uplo, f77char *transa, f77char *diag, f77integer *m, f77integer *n, f77doublecomplex *alpha, f77doublecomplex *a, f77integer *lda, f77doublecomplex *b, f77integer *ldb) = &ZTRSM;
f77integer (* dtrsm)(f77char *side, f77char *uplo, f77char *transa, f77char *diag, f77integer *m, f77integer *n, f77doublereal *alpha, f77doublereal *a, f77integer *lda, f77doublereal *b, f77integer *ldb) = &DTRSM;
f77integer (* sswap)(f77integer *n, f77real *sx, f77integer *incx, f77real *sy, f77integer *incy) = &SSWAP;
f77integer (* zcopy)(f77integer *n, f77doublecomplex *zx, f77integer *incx, f77doublecomplex *zy, f77integer *incy) = &ZCOPY;
f77integer (* dtbmv)(f77char *uplo, f77char *trans, f77char *diag, f77integer *n, f77integer *k, f77doublereal *a, f77integer *lda, f77doublereal *x, f77integer *incx) = &DTBMV;
f77integer (* chpr)(f77char *uplo, f77integer *n, f77real *alpha, f77complex *x, f77integer *incx, f77complex *ap) = &CHPR;
f77integer (* dspr)(f77char *uplo, f77integer *n, f77doublereal *alpha, f77doublereal *x, f77integer *incx, f77doublereal *ap) = &DSPR;
f77doublereal (* dasum)(f77integer *n, f77doublereal *dx, f77integer *incx) = &DASUM;
f77doublereal (* dcabs1)(f77doublecomplex *z__) = &DCABS1;
f77integer (* sger)(f77integer *m, f77integer *n, f77real *alpha, f77real *x, f77integer *incx, f77real *y, f77integer *incy, f77real *a, f77integer *lda) = &SGER;
f77integer (* dcopy)(f77integer *n, f77doublereal *dx, f77integer *incx, f77doublereal *dy, f77integer *incy) = &DCOPY;
};