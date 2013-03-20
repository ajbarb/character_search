#define MAXN 2048             /* max 2d fft size */
#define EPSILON 0.00001       /* for comparing fp numbers */
#define PI 3.14159265358979   /* 4*atan(1.0) */

typedef struct {float r,i;} mycomplex; 

/* swap a pair of complex numbers */
#define SWAP(a,b) {float swap_temp=(a).r;(a).r=(b).r;(b).r=swap_temp;\
		         swap_temp=(a).i;(a).i=(b).i;(b).i=swap_temp;}
