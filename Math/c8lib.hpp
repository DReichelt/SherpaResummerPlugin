double            c8_abs ( Complex x );
Complex  c8_acos ( Complex c1 );
Complex  c8_acosh ( Complex c1 );
Complex  c8_add ( Complex c1, Complex c2 );
double            c8_arg ( Complex x );
Complex  c8_asin ( Complex c1 );
Complex  c8_asinh ( Complex c1 );
Complex  c8_atan ( Complex c1 );
Complex  c8_atanh ( Complex c1 );
Complex  c8_conj ( Complex c1 );
void              c8_copy ( Complex c1, Complex c2 );
Complex  c8_cos ( Complex c1 );
Complex  c8_cosh ( Complex c1 );
Complex  c8_cube_root ( Complex x );
Complex  c8_div ( Complex c1, Complex c2 );
Complex  c8_div_r8 ( Complex c1, double r );
Complex  c8_exp ( Complex c1 );
Complex  c8_i ( );
double            c8_imag ( Complex c );
Complex  c8_inv ( Complex c1 );
bool              c8_le_l1 ( Complex x, Complex y );
bool              c8_le_l2 ( Complex x, Complex y );
bool              c8_le_li ( Complex x, Complex y );
Complex  c8_log ( Complex c1 );
double            c8_mag ( Complex x );
Complex  c8_mul ( Complex c1, Complex c2 );
Complex  c8_neg ( Complex c1 );
Complex  c8_nint ( Complex c1 );
double            c8_norm_l1 ( Complex x );
double            c8_norm_l2 ( Complex x );
double            c8_norm_li ( Complex x );
Complex  c8_normal_01 ( int *seed );
Complex  c8_one ( );
void              c8_print ( Complex a, string title );
double            c8_real ( Complex c );
Complex  c8_sin ( Complex c1 );
Complex  c8_sinh ( Complex c1 );
Complex  c8_sqrt ( Complex x );
Complex  c8_sub ( Complex c1, Complex c2 );
void              c8_swap ( Complex *x, Complex *y );
Complex  c8_tan ( Complex c1 );
Complex  c8_tanh ( Complex c1 );
void              c8_to_cartesian ( Complex c, double *x, double *y );
void              c8_to_polar ( Complex c, double *r, double *theta );
Complex  c8_uniform_01 ( int *seed );
Complex  c8_zero ( );
void              c8mat_add ( int m, int n, Complex alpha, Complex a[],
                  Complex beta, Complex b[], Complex c[] );
void              c8mat_add_r8 ( int m, int n, double alpha, Complex a[],
                  double beta, Complex b[], Complex c[] );
void              c8mat_copy ( int m, int n, Complex a1[], 
                  Complex a2[] );
Complex *c8mat_copy_new ( int m, int n, Complex a1[] );
void              c8mat_fss ( int n, Complex a[], int nb, Complex x[] );
Complex *c8mat_fss_new ( int n, Complex a[], int nb, 
                  Complex b[] );
Complex *c8mat_identity_new ( int n );
Complex *c8mat_indicator_new ( int m, int n );
void              c8mat_minvm ( int n1, int n2, Complex a[], 
                  Complex b[], Complex e[] );
Complex *c8mat_minvm_new ( int n1, int n2, Complex a[], 
                  Complex b[] );
void              c8mat_mm ( int n1, int n2, int n3, Complex a[], 
                  Complex b[], Complex c[] );
Complex *c8mat_mm_new ( int n1, int n2, int n3, Complex a[], 
                  Complex b[] );
void              c8mat_nint ( int m, int n, Complex a[] );
double            c8mat_norm_fro ( int m, int n, Complex a[] );
double            c8mat_norm_l1 ( int m, int n, Complex a[] );
double            c8mat_norm_li ( int m, int n, Complex a[] );
void              c8mat_print ( int m, int n, Complex a[], string title );
void              c8mat_print_some ( int m, int n, Complex a[], int ilo, 
                  int jlo, int ihi, int jhi, string title );
void              c8mat_scale ( int m, int n, Complex alpha, Complex a[] );
void              c8mat_scale_r8 ( int m, int n, double alpha, Complex a[] );
void              c8mat_uniform_01 ( int m, int n, int *seed, 
                  Complex c[] );
Complex *c8mat_uniform_01_new ( int m, int n, int *seed );
Complex *c8mat_zero_new ( int m, int n );
void              c8vec_copy ( int n, Complex a1[], 
                  Complex a2[] );
Complex *c8vec_copy_new ( int n, Complex a1[] );
Complex *c8vec_indicator_new ( int n );
void              c8vec_nint ( int n, Complex a[] );
double            c8vec_norm_l2 ( int n, Complex a[] );
void              c8vec_print ( int n, Complex a[], string title );
void              c8vec_print_part ( int n, Complex a[], int max_print, 
                  string title );
void              c8vec_print_some ( int n, Complex a[], int i_lo, 
                  int i_hi, string title );
void              c8vec_sort_a_l2 ( int n, Complex x[] );
Complex *c8vec_spiral ( int n, int m, Complex c1, 
                  Complex c2 );
Complex *c8vec_uniform_01_new ( int n, int *seed );
Complex *c8vec_unity ( int n );
Complex  cartesian_to_c8 ( double x, double y );
Complex  polar_to_c8 ( double r, double theta );
Complex  r8_csqrt ( double x );
double            r8_uniform_01 ( int *seed );
void              r8poly2_root ( double a, double b, double c, 
                  Complex *r1, Complex *r2 );
void              r8poly3_root ( double a, double b, double c, double d, 
                  Complex *r1, Complex *r2, 
                  Complex *r3 );
void              r8poly4_root ( double a, double b, double c, double d, double e,
                  Complex *r1, Complex *r2, Complex *r3,
                  Complex *r4 );
void              sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
