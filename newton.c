#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <math.h>

void complex_power();
void calc_roots();
void newton();

int main(int argc, char** argv) {
  if (argc != 4) {
    printf("Wrong number of arguments\n");
    return 0;
  }
  int c;
  int tval;
  int lval;
  int exp;
  while ((c = getopt(argc, argv, "t:l:"))!= -1) {
    switch (c) {
    case 't':
      tval = atoi(optarg);
    case 'l':
      lval = atoi(optarg);
    }
  }
  exp = atoi(argv[3]);
  
  double x0_re = -3.0;
  double x0_im = -1.0;
  int root;
  double roots_re[exp];
  double roots_im[exp];
  calc_roots(roots_re, roots_im, &exp);
  newton(&x0_re, &x0_im, &root, &roots_re, &roots_im, &exp);

  if (root <= exp) {
    printf("With x_0 = %.2f + %.2fi the root %.2f + %.2fi was found\n", x0_re, x0_im, roots_re[root], roots_im[root]);
  }
  else {
    printf("Newton's method did not converge\n");
  }
  
}


void calc_roots(double* roots_re, double* roots_im, int* exp) {

  for (size_t i = 0; i < *exp; i++) {
    roots_re[i] = cos(2*M_PI*i / *exp);
    roots_im[i] = sin(2*M_PI*i / *exp);
    
  }
  
}

void newton(double* x0_re, double* x0_im, int* root, double* roots_re, double* roots_im, int* exp) {
  #define MAX 10000000000
  #define MIN 0.000001
  //double roots_re[*exp];
  //double roots_im[*exp];
  //calc_roots(roots_re, roots_im, exp); 
  double xk_re = *x0_re;
  double xk_im = *x0_im;

  /* Perform newton */
  int found_root = 0;
  int divergent = 0;
  int origin = 0;

  /* Next iteration value */
  double xk1_re;
  double xk1_im;

  /* Function value in curr point */
  double fxk_re;
  double fxk_im;

  /* Function derivative value in curr point */
  double fpxk_re;
  double fpxk_im;

  /* Exponent in derivative */
  int pexp = *exp-1;
  for (size_t k = 0; k < 100; k++) {
    divergent = xk_re>MAX || -xk_re>MAX || xk_im>MAX || -xk_im>MAX;
    origin = xk_re*xk_re + xk_im*xk_im < MIN;
    if (divergent || origin) {
     *root = *exp + 1;
     break;
    }
    if (found_root) {
      break;
    }
    //printf("root not found, cont with x_%i = %f + %fi\n", k, xk_re, xk_im);
    /*Check if close enough to real root*/
    double dist;
    for (size_t i = 0; i < *exp; i++) {
      dist = (xk_re-roots_re[i])*(xk_re-roots_re[i]) + (xk_im-roots_im[i])*(xk_im-roots_im[i]);
      if (dist < MIN) {
	printf("found root! xk_re = %f xk_im = %f\n", xk_re, xk_im);
	*root= i;
	found_root = 1;
      }
    }
    
    complex_power(&fxk_re, &fxk_im, &xk_re, &xk_im, exp);
    
    fxk_re = fxk_re -1;
    complex_power(&fpxk_re, &fpxk_im, &xk_re, &xk_im, &pexp);
    fpxk_re *= *exp;
    fpxk_im *= *exp;
    xk1_re = xk_re - ((fxk_re*fpxk_re + fxk_im*fpxk_im) / (fpxk_re*fpxk_re + fpxk_im*fpxk_im));
    xk1_im = xk_im - (fxk_im*fpxk_re - fxk_re*fpxk_im) / (fpxk_re*fpxk_re + fpxk_im*fpxk_im);
    /* Update current root */
    xk_re = xk1_re;
    xk_im = xk1_im;
  }
}




/* Computes the power of a complex number */
void complex_power(double* res_re, double* res_im, double* x_re, double* x_im, int* exp) {
  int i = 0;
  double tmp_re;
  double tmp_im;
  *res_re = *x_re;
  *res_im = *x_im;
  tmp_re = *x_re;
  tmp_im = *x_im;
  
  while (i < *exp-1) {
    *res_re = (tmp_re) * (*x_re) - (tmp_im) * (*x_im);
    *res_im = (tmp_im) * (*x_re) + (tmp_re) * (*x_im);
    tmp_re = *res_re;
    tmp_im = *res_im;
    i++;
  }
}
