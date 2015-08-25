#include <stdio.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_odeiv.h>

#define INITIAL_R 1e-6
#define RK_ACCURACY 1e-6
#define BIG_REAL 1e6

static int dim;
static int rank;
static double *s, *m, *ss, *mm, *x, *y;
static double radius=20.0; /* radius=20.0 or 1000.0 */
static double Log_C;

static void ball_prob(void);
static void ball_prob2(void);
static void ball_prob_graph(void);
static void fisher_bingham(void);
static void fisher_bingham_graph(void);
static void hirotsu_1_1(void);
static void hirotsu_1_2(void);
static void hirotsu_2_1(void);
static void hirotsu_2_2(void);
static void anderson_darling(void);
static double gamma_function1(int);
static double gamma_function2(int);
static void comparison(void);
static void feller(void);

static int function_1(double r, const double f[],  double df[], void *prms);
static int function_2(double r, const double f[],  double df[], void *prms);
static int function_3(double r, const double f[],  double df[], void *prms);
static int function_4(double r, const double f[],  double df[], void *prms);
static void runge_kutta_1(double *f, double r0, double r1);
static void runge_kutta_2(double *f, double r0, double r1);
static void runge_kutta_3(double *f, double r0, double r1);
static void runge_kutta_4(double *f, double r0, double r1);
static void show_ratio(double r, double *f);
static void print_vector(FILE *fp, double *vector, int length, char *string);

static void 
ball_prob(void)
{
  /* Set parameters */
  rank = 2 * dim + 1;
  int i;
  for ( i = 0; i < dim; i++){
    ss[i] = s[i] * s[i];
    mm[i] = m[i] * m[i];
    x[i] = -0.5 / (s[i]*s[i]);
    y[i] = m[i] / (s[i]*s[i]);    
  }
  /*
    print_vector(stdout, s, dim, "  sigma:");
    print_vector(stdout, m, dim, "     mu:");
    print_vector(stdout, x, dim, "      x:");
    print_vector(stdout, y, dim, "      y:");
  */

  /* initial value */
  double f[rank];
  for (i = 0; i < dim; i++){
    f[i]     = y[i];
    f[i+dim] = 1.0;
  }
  f[rank-1] = 0.0;
  //  print_vector(stdout, f, rank, "      f:");
  double r = INITIAL_R;

  /* set Log_C*/
  Log_C= (dim+1)*log(r);
  for (i = 0; i < dim; i++)    Log_C -= log(s[i]);
  for (i = 0; i < dim; i++)    Log_C -= 0.5*mm[i]/ss[i];
  for (i = dim; i > 0; i -= 2) Log_C -= log(i);
  if (dim % 2) Log_C += 0.5 * (log(2.0) - log(M_PI));
  //  printf("Log_C: %.10g\n", Log_C);
  //  printf("    C: %.10g\n", exp(Log_C));

  if (radius < 1.0) {
    /* move r form 1e-6 to radius */  
    runge_kutta_1(f, INITIAL_R, radius);
    print_vector(stdout, f, rank, "f:\t");
    r = radius;
  } else {
    /* move r form 1e-6 to 1.0 */  
    runge_kutta_1(f, INITIAL_R, 1.0);
    r = 1.0;
    
    //  fprintf(stdout, "r=%f, C=%g\n", r, exp(Log_C));
    //  print_vector(stdout, f, rank, "f:\t");
    
    for ( i = 0 ; i < rank-1; i++)
      f[i] *= exp(-r*r*x[0]-r*y[0]);
    f[0] /= r;
    f[dim] /= r*r;
    //print_vector(stdout, f, rank, "f:\t");

    /* move r form 1.0 to radius */
    runge_kutta_2(f, 1.0, radius);
    r = radius;
  }

  //fprintf(stdout, "r=%f\n", r);
  //print_vector(stdout, f, rank, "f:\t");

  fprintf(stdout, "Probability= %10.6f\n", f[rank-1]);

  return;
}

static void 
ball_prob2(void)
{
  /* Set parameters */
  rank = 2 * dim + 1;
  int i;
  for ( i = 0; i < dim; i++){
    ss[i] = s[i] * s[i];
    mm[i] = m[i] * m[i];
    x[i] = -0.5 / (s[i]*s[i]);
    y[i] = m[i] / (s[i]*s[i]);    
  }
  /*
    print_vector(stdout, s, dim, "  sigma:");
    print_vector(stdout, m, dim, "     mu:");
    print_vector(stdout, x, dim, "      x:");
    print_vector(stdout, y, dim, "      y:");
  */

  /* initial value */
  double f[rank];
  for (i = 0; i < dim; i++){
    f[i]     = y[i];
    f[i+dim] = 1.0;
  }
  f[rank-1] = 0.0;
  //  print_vector(stdout, f, rank, "      f:");
  double r = INITIAL_R;

  /* set Log_C*/
  Log_C= (dim+1)*log(r);
  for (i = 0; i < dim; i++)    Log_C -= log(s[i]);
  for (i = 0; i < dim; i++)    Log_C -= 0.5*mm[i]/ss[i];
  for (i = dim; i > 0; i -= 2) Log_C -= log(i);
  if (dim % 2) Log_C += 0.5 * (log(2.0) - log(M_PI));
  //  printf("Log_C: %.10g\n", Log_C);
  //  printf("    C: %.10g\n", exp(Log_C));

  if (radius < 1.0) {
    /* move r form 1e-6 to radius */  
    runge_kutta_1(f, INITIAL_R, radius);
    print_vector(stdout, f, rank, "f:\t");
    r = radius;
  } else {
    /* move r form 1e-6 to 1.0 */  
    runge_kutta_1(f, INITIAL_R, 1.0);
    r = 1.0;
    
    //  fprintf(stdout, "r=%f, C=%g\n", r, exp(Log_C));
    //  print_vector(stdout, f, rank, "f:\t");
    
    for ( i = 0 ; i < rank-1; i++)
      f[i] *= exp(-r*r*x[0]-r*y[0]);
    f[0] /= r;
    f[dim] /= r*r;
    //print_vector(stdout, f, rank, "f:\t");

    /* move r form 1.0 to radius */
    runge_kutta_2(f, 1.0, radius);
    r = radius;
  }

  //fprintf(stdout, "r=%f\n", r);
  //print_vector(stdout, f, rank, "f:\t");

  //fprintf(stdout, "Probability= %10.6f\n", f[rank-1]);
  fprintf(stdout, "p-1= %4.2e\n", f[rank-1]-1);

  return;
}

static void 
ball_prob_graph(void)
{
  /* Set parameters */
  rank = 2 * dim + 1;
  int i;
  for ( i = 0; i < dim; i++){
    ss[i] = s[i] * s[i];
    mm[i] = m[i] * m[i];
    x[i] = -0.5 / (s[i]*s[i]);
    y[i] = m[i] / (s[i]*s[i]);    
  }

  /* initial value */
  double f[rank];
  for (i = 0; i < dim; i++){
    f[i]     = y[i];
    f[i+dim] = 1.0;
  }
  f[rank-1] = 0.0;
  double r = INITIAL_R;

  /* set Log_C*/
  Log_C= (dim+1)*log(r);
  for (i = 0; i < dim; i++)    Log_C -= log(s[i]);
  for (i = 0; i < dim; i++)    Log_C -= 0.5*mm[i]/ss[i];
  for (i = dim; i > 0; i -= 2) Log_C -= log(i);
  if (dim % 2) Log_C += 0.5 * (log(2.0) - log(M_PI));

  int k = 0;
  double from, to;
  from = INITIAL_R;
  to = (k+1) * 0.01 * radius;
  while ( to < 1.0 ){
    runge_kutta_1(f, from, to);
    fprintf(stdout, "%g \t %g\n", to, f[rank-1]);
    k++;
    from = to;
    to = (k+1) * 0.01 * radius;
  }

  runge_kutta_1(f, from, 1.0);
  r = 1.0;

  for ( i = 0 ; i < rank-1; i++)
    f[i] *= exp(-r*r*x[0]-r*y[0]);
  f[0] /= r;
  f[dim] /= r*r;

  from = 1.0;

  while ( k < 100 ){
    runge_kutta_2(f, from, to);
    fprintf(stdout, "%g \t %g\n", to, f[rank-1]);
    k++;
    from = to;
    to = (k+1) * 0.01 * radius;
  }

  return;
}

static void 
fisher_bingham(void)
{
  /* Set parameters */
  rank = 2 * dim;
  int i;
  for ( i = 0; i < dim; i++){
    ss[i] = s[i] * s[i];
    mm[i] = m[i] * m[i];
    x[i] = -0.5 / (s[i]*s[i]);
    y[i] = m[i] / (s[i]*s[i]);    
  }
  /*
  print_vector(stdout, s, dim, "  sigma:");
  print_vector(stdout, m, dim, "     mu:");
  print_vector(stdout, x, dim, "      x:");
  print_vector(stdout, y, dim, "      y:");
  */

  /* initial value */
  double f[rank];
  for (i = 0; i < dim; i++){
    f[i]     = y[i];
    f[i+dim] = 1.0;
  }
  //f[rank-1] = 0.0;
  //print_vector(stdout, f, rank, "      f:");
  double r = INITIAL_R;

  /* set Log_C*/
  Log_C= (dim+1)*log(r);
  for (i = dim; i > 0; i -= 2) Log_C -= log(i);
  if (dim % 2) {
    Log_C += log(2.0) + 0.5*(dim-1)*log(2*M_PI);
  } else {
    Log_C += 0.5*dim*log(2*M_PI);
  }
  //  printf("Log_C: %.10g\n", Log_C);
  //  printf("    C: %.10g\n", exp(Log_C));

  /*
  double max_x = x[0];
  printf("  x_1: %.10g\n", max_x);
  printf("   -r^2*x_1 : %.10g\n",     -r*r*max_x);
  printf("e^(-r^2*x_1): %.10g\n", exp(-r*r*max_x));
  */

  /* move r form 1e-6 to 1.0 */  
  runge_kutta_3(f, INITIAL_R, 1.0);
  r = 1.0;

  //fprintf(stdout, "r=%f, C=%g\n", r, exp(Log_C));
  //print_vector(stdout, f, rank, "f:\t");

  for ( i = 0 ; i < rank; i++)
    f[i] *= exp(-r*r*x[0]-r*y[0]);
  f[0] /= r;
  f[dim] /= r*r;
  //print_vector(stdout, f, rank, "f:\t");

  /* move r form 1.0 to radius */
  runge_kutta_4(f, 1.0, radius);
  r = radius;

  //fprintf(stdout, "r=%f\n", r);
  //print_vector(stdout, f, rank, "f:\t");

  /* Laprace approximations */
  double large_num1, large_num2;
  large_num1= r*r*x[0] + log( exp(r*y[0])+exp(-r*y[0]) );
  large_num2=Log_C+r*r*x[0]+r*y[0];

  double sum, prod;
  sum = 0.0;
  for ( i = 1 ; i < dim; i++)  
    sum -= y[i]*y[i] / ( 4*(x[i]-x[0]) );
  prod = pow(M_PI, (dim-1)*0.5) * exp(sum);
  for ( i = 1 ; i < dim; i++)  
    prod *= 1.0 / pow( x[0]-x[i], 0.5);
  double lap[rank];
  lap[0] = r*(exp(r*y[0])-exp(-r*y[0]))/(exp(r*y[0]) + exp(-r*y[0]))*prod;
  lap[dim] = r*r*prod;
  for ( i = 1 ; i < dim; i++){
    double c = 0.5 / (x[0]-x[i]);
    lap[i] = y[i]*c*prod;
    lap[i+dim] = (y[i]*y[i]*c*c+c)*prod;
  }
  //print_vector(stdout, lap, rank, "lap:\t");
  //fprintf(stdout, "\t times exp(%lf)\n", large_num1);

  double g[rank];
  /* Differantials of FB integral */
  for ( i = 0 ; i < rank; i++)
    g[i] = f[i] * exp(large_num2);
  g[0] *= r;
  g[dim] *= r*r;
  print_vector(stdout, g, rank, "f:\t");

  /* laprace approximations */
  for ( i = 0 ; i < rank; i++)
    g[i] = lap[i] * exp(large_num1);
  print_vector(stdout, g, rank, "lap:\t");

  /* (Differantial of FB integral) / (Laprace approximation) */
  for ( i = 0 ; i < rank; i++)
    g[i] = f[i]/lap[i] * exp(large_num2-large_num1);
  g[0] *= r;
  g[dim] *= r*r;
  print_vector(stdout, g, rank, "ratio:\t");

  /* (Fisher-Bingham integral) / (Laprace apporoximation) */
  double fb;
  fb = r*r*f[dim];
  for ( i = 1; i < dim; i++)
    fb += f[dim+i];
  fb /= r*r;
  fprintf(stdout, "Fisher-Bingham integral:\t%g\n", fb*exp(large_num2) );
  fprintf(stdout, "  Laprace approximation:\t%g\n", prod*exp(large_num1));
  
  fprintf(stdout,"(FB)/(Laprace):\t\t\t%f\n",fb/prod*exp(large_num2-large_num1));

  return;
}

static void 
fisher_bingham_graph(void)
{ 
  /* Set parameters */
  rank = 2 * dim;
  int i;
  for ( i = 0; i < dim; i++){
    ss[i] = s[i] * s[i];
    mm[i] = m[i] * m[i];
    x[i] = -0.5 / (s[i]*s[i]);
    y[i] = m[i] / (s[i]*s[i]);    
  }

  /* initial value */
  double f[rank];
  for (i = 0; i < dim; i++){
    f[i]     = y[i];
    f[i+dim] = 1.0;
  }
  double r = INITIAL_R;

  /* set Log_C*/
  Log_C= (dim+1)*log(r);
  for (i = dim; i > 0; i -= 2) Log_C -= log(i);
  if (dim % 2) {
    Log_C += log(2.0) + 0.5*(dim-1)*log(2*M_PI);
  } else {
    Log_C += 0.5*dim*log(2*M_PI);
  }

  int k = 0;
  int n_points = 100;
  double step_size = radius / n_points;
  double from, to;
  from = INITIAL_R;
  to = (k+1) * step_size;
  while ( to < 1.0 ){
    runge_kutta_3(f, from, to);
    show_ratio(to, f);
    k++;
    from = to;
    to = (k+1) * step_size;
  }

  runge_kutta_3(f, from, 1.0);

  r = 1.0;
  for ( i = 0 ; i < rank; i++)
    f[i] *= exp(-r*r*x[0]-r*y[0]);
  f[0] /= r;
  f[dim] /= r*r;

  from = 1.0;

  while ( k < n_points ){
    runge_kutta_4(f, from, to);
    show_ratio(to, f);
    k++;
    from = to;
    to = (k+1) * step_size;
  }
  return;
}

static void
hirotsu_1_1(void)
{
  double sigma_mu[6*dim];
  s  = sigma_mu;
  m  = sigma_mu +   dim;
  ss = sigma_mu + 2*dim;
  mm = sigma_mu + 3*dim;
  x  = sigma_mu + 4*dim;
  y  = sigma_mu + 5*dim;

  int i;
  for (i=0; i<dim; i++){
    s[i] = sqrt( (double) (dim+1) / ((i+1)*(i+2)) );
    m[i] = 0.0;
  }

  ball_prob2();
  return;
}

static void
hirotsu_1_2(void)
{
  double sigma_mu[6*dim];
  s  = sigma_mu;
  m  = sigma_mu +   dim;
  ss = sigma_mu + 2*dim;
  mm = sigma_mu + 3*dim;
  x  = sigma_mu + 4*dim;
  y  = sigma_mu + 5*dim;

  int i;
  for (i=0; i<dim; i++){
    s[i] = sqrt( (double) (dim+1) / ((i+1)*(i+2)) );
    m[i] = 0.01 * i;
  }

  ball_prob2();
  //fisher_bingham();
  return;
}

static void
hirotsu_2_1(void)
{
  double sigma_mu[6*dim];
  s  = sigma_mu;
  m  = sigma_mu +   dim;
  ss = sigma_mu + 2*dim;
  mm = sigma_mu + 3*dim;
  x  = sigma_mu + 4*dim;
  y  = sigma_mu + 5*dim;

  int i;
  for (i=0; i<dim; i++){
    s[i] = sqrt( (double) 2.0 * (dim+2) * (dim+3) / ((i+1)*(i+2)*(i+3)*(i+4)) );
    m[i] = 0.0;
  }

  ball_prob2();
  return;
}

static void
hirotsu_2_2(void)
{
  double sigma_mu[6*dim];
  s  = sigma_mu;
  m  = sigma_mu +   dim;
  ss = sigma_mu + 2*dim;
  mm = sigma_mu + 3*dim;
  x  = sigma_mu + 4*dim;
  y  = sigma_mu + 5*dim;

  int i;
  for (i=0; i<dim; i++){
    s[i] = sqrt( (double) 2.0 * (dim+2) * (dim+3) / ((i+1)*(i+2)*(i+3)*(i+4)) );
    m[i] = 0.01 * i;
  }

  ball_prob2();
  return;
}

static void 
anderson_darling(void)
{
  double sigma_mu[6*dim];
  s  = sigma_mu;
  m  = sigma_mu +   dim;
  ss = sigma_mu + 2*dim;
  mm = sigma_mu + 3*dim;
  x  = sigma_mu + 4*dim;
  y  = sigma_mu + 5*dim;

  int i;
  for (i=0; i<dim; i++){
    s[i] = sqrt( (double) 1.0 / ((i+1)*(i+2)) );
    m[i] = 0.0;
  }

  ball_prob2();
  return;
}

static int 
function_1(double r, const double f[],  double df[], void *prms)
{
  const double  *f2 = f  + dim;
  double sum = 0.0;
  double inv_r = 1.0 / r;
  int i;
  for (i = 0; i < dim; i++)    
    sum += f2[i];
  for (i = 0; i < dim; i++)
    df[i]  = 2*r*x[i]*f[i] + inv_r * (f[i]+y[i]*sum);
  for (i = 0; i < dim; i++){
    df[i+dim] = r*(y[i]*f[i] + 2*x[i]*f2[i]) ;
    df[i+dim] += inv_r*(f2[i]+sum);
  }

  df[rank-1] = inv_r * inv_r * sum * exp(Log_C);

  return GSL_SUCCESS;
}


static int 
function_2(double r, const double f[],  double df[], void *prms)
{
  const double  *f2 = f  + dim;
  double g[2*dim];
  double *g2 = g  + dim;
  double sum = 0.0;
  double inv_r = 1.0 / r, rr = r*r;
  int i;

  g[0] = r * f[0];
  for (i = 1; i < dim; i++)    
    g[i] = f[i];
  g2[0] = rr * f2[0];
  for (i = 1; i < dim; i++)    
    g2[i] = f2[i];

  for (i = 0; i < dim; i++)    
    sum += g2[i];
  for (i = 0; i < dim; i++)
    df[i]  = 2*r*x[i]*g[i] + inv_r * (g[i]+y[i]*sum);
  for (i = 0; i < dim; i++){
    df[i+dim] = r*(y[i]*g[i] + 2*x[i]*g2[i]) ;
    df[i+dim] += inv_r*(g2[i]+sum);
  }

  df[0] *= inv_r;
  df[dim] *= 1.0/rr;

  double m2_rx1 = -2.0*r*x[0]-y[0]; 
  for (i = 0; i < 2*dim; i++)
    df[i]  += m2_rx1 * f[i];
  df[0] += -inv_r*f[0];
  df[dim] += -2.0*inv_r*f[dim];

  df[rank-1] = exp(Log_C+rr*x[0]+r*y[0]) * sum  / rr;

  return GSL_SUCCESS;
}

static int 
function_3(double r, const double f[],  double df[], void *prms)
{
  const double  *f2 = f  + dim;
  double sum = 0.0;
  double inv_r = 1.0 / r;
  int i;
  for (i = 0; i < dim; i++)    
    sum += f2[i];
  for (i = 0; i < dim; i++)
    df[i]  = 2*r*x[i]*f[i] + inv_r * (f[i]+y[i]*sum);
  for (i = 0; i < dim; i++){
    df[i+dim] = r*(y[i]*f[i] + 2*x[i]*f2[i]) ;
    df[i+dim] += inv_r*(f2[i]+sum);
  }

  return GSL_SUCCESS;
}


static int 
function_4(double r, const double f[],  double df[], void *prms)
{
  const double  *f2 = f  + dim;
  double g[2*dim];
  double *g2 = g  + dim;
  double sum = 0.0;
  double inv_r = 1.0 / r, rr = r*r;
  int i;

  g[0] = r * f[0];
  for (i = 1; i < dim; i++)    
    g[i] = f[i];
  g2[0] = rr * f2[0];
  for (i = 1; i < dim; i++)    
    g2[i] = f2[i];

  for (i = 0; i < dim; i++)    
    sum += g2[i];
  for (i = 0; i < dim; i++)
    df[i]  = 2*r*x[i]*g[i] + inv_r * (g[i]+y[i]*sum);
  for (i = 0; i < dim; i++){
    df[i+dim] = r*(y[i]*g[i] + 2*x[i]*g2[i]) ;
    df[i+dim] += inv_r*(g2[i]+sum);
  }

  df[0] *= inv_r;
  df[dim] *= 1.0/rr;

  double m2_rx1 = -2.0*r*x[0]-y[0]; 
  for (i = 0; i < 2*dim; i++)
    df[i]  += m2_rx1 * f[i];
  df[0] += -inv_r*f[0];
  df[dim] += -2.0*inv_r*f[dim];

  return GSL_SUCCESS;
}

static void 
runge_kutta_1(double *f, double r0, double r1)
{
  double h = RK_ACCURACY;
  const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
  gsl_odeiv_step * s  = gsl_odeiv_step_alloc (T, rank);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new (h, 0.0);
  gsl_odeiv_evolve * e  = gsl_odeiv_evolve_alloc (rank);
  gsl_odeiv_system sys1 = {function_1, NULL, rank, NULL};

  double r = r0;
  double cnst = 10, inv_cnst = 1.0/cnst, log_cnst = log(cnst);
  while (r < r1){
    int status = gsl_odeiv_evolve_apply (e, c, s, &sys1, &r, r1, &h, f);
    if (status != GSL_SUCCESS)
      break;
    if (f[dim] > BIG_REAL){ /* <- changed! rank-1 --> dim */
      int i;
      for (i = 0; i < rank-1; i++) f[i] *= inv_cnst;
      Log_C += log_cnst;
    }
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);

  return;
}

static void 
runge_kutta_2(double *f, double r0, double r1)
{
  double h = RK_ACCURACY;
  const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
  gsl_odeiv_step * s  = gsl_odeiv_step_alloc (T, rank);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new (h, 0.0);
  gsl_odeiv_evolve * e  = gsl_odeiv_evolve_alloc (rank);
  gsl_odeiv_system sys2 = {function_2, NULL, rank, NULL};

  double r = r0;
  double cnst = 10, inv_cnst = 1.0/cnst, log_cnst = log(cnst);
  while (r < r1){
    int status = gsl_odeiv_evolve_apply (e, c, s, &sys2, &r, r1, &h, f);
    if (status != GSL_SUCCESS)
      break;
    if (f[dim] > BIG_REAL){ /* <- changed! rank-1 --> dim */
      int i;
      for (i = 0; i < rank-1; i++) f[i] *= inv_cnst;
      Log_C += log_cnst;
    }
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);

  return;
}

static void 
runge_kutta_3(double *f, double r0, double r1)
{
  double h = RK_ACCURACY;
  const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
  gsl_odeiv_step * s  = gsl_odeiv_step_alloc (T, rank);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new (h, 0.0);
  gsl_odeiv_evolve * e  = gsl_odeiv_evolve_alloc (rank);
  gsl_odeiv_system sys3 = {function_3, NULL, rank, NULL};

  double r = r0;
  double cnst = 10, inv_cnst = 1.0/cnst, log_cnst = log(cnst);
  while (r < r1){
    int status = gsl_odeiv_evolve_apply (e, c, s, &sys3, &r, r1, &h, f);
    if (status != GSL_SUCCESS)
      break;
    if (f[dim] > BIG_REAL){ /* <- changed! rank-1 --> dim */
      int i;
      for (i = 0; i < 2*dim; i++) f[i] *= inv_cnst;
      Log_C += log_cnst;
    }
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);

  return;
}

static void 
runge_kutta_4(double *f, double r0, double r1)
{
  double h = RK_ACCURACY;
  const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
  gsl_odeiv_step * s  = gsl_odeiv_step_alloc (T, rank);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new (h, 0.0);
  gsl_odeiv_evolve * e  = gsl_odeiv_evolve_alloc (rank);
  gsl_odeiv_system sys4 = {function_4, NULL, rank, NULL};

  double r = r0;
  double cnst = 10, inv_cnst = 1.0/cnst, log_cnst = log(cnst);
  while (r < r1){
    int status = gsl_odeiv_evolve_apply (e, c, s, &sys4, &r, r1, &h, f);
    if (status != GSL_SUCCESS)
      break;
    if (f[dim] > BIG_REAL){ /* <- changed! rank-1 --> dim */
      int i;
      for (i = 0; i < rank; i++) f[i] *= inv_cnst;
      Log_C += log_cnst;
    }
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);

  return;
}

static void 
show_ratio(double r, double *f)
{
  fprintf(stdout, "%10.6g ", r);

  /* Laprace approximations */
  double large_num1, large_num2;
  large_num1= r*r*x[0] + log( exp(r*y[0])+exp(-r*y[0]) );
  large_num2=Log_C+r*r*x[0]+r*y[0];

  int i;
  double sum, prod;
  sum = 0.0;
  for ( i = 1 ; i < dim; i++)  
    sum -= y[i]*y[i] / ( 4*(x[i]-x[0]) );
  prod = pow(M_PI, (dim-1)*0.5) * exp(sum);
  for ( i = 1 ; i < dim; i++)  
    prod *= 1.0 / pow( x[0]-x[i], 0.5);
  double lap[rank];
  lap[0] = r*(exp(r*y[0])-exp(-r*y[0]))/(exp(r*y[0]) + exp(-r*y[0]))*prod;
  lap[dim] = r*r*prod;
  for ( i = 1 ; i < dim; i++){
    double c = 0.5 / (x[0]-x[i]);
    lap[i] = y[i]*c*prod;
    lap[i+dim] = (y[i]*y[i]*c*c+c)*prod;
  }
  //print_vector(stdout, lap, rank, "lap:\t");
  //fprintf(stdout, "\t times exp(%lf)\n", large_num1);

  double g[rank];
  /* Differantials of FB integral */
  /*
  for ( i = 0 ; i < rank; i++)
    g[i] = f[i] * exp(large_num2);
  g[0] *= r;
  g[dim] *= r*r;
  print_vector(stdout, g, rank, "f:\t");
  */

  /* laprace approximations */
  /*
  for ( i = 0 ; i < rank; i++)
    g[i] = lap[i] * exp(large_num1);
  print_vector(stdout, g, rank, "lap:\t");
  */

  /* (Differantial of FB integral) / (Laprace approximation) */
  for ( i = 0 ; i < rank; i++)
    g[i] = f[i]/lap[i] * exp(large_num2-large_num1);
  g[0] *= r;
  g[dim] *= r*r;
  //print_vector(stdout, g, rank, "ratio:\t");

  /* (Fisher-Bingham integral) / (Laprace apporoximation) */
  double fb;
  fb = r*r*f[dim];
  for ( i = 1; i < dim; i++)
    fb += f[dim+i];
  fb /= r*r;
  //fprintf(stdout, "Fisher-Bingham integral:\t%g\n", fb*exp(large_num2) );
  //fprintf(stdout, "  Laprace approximation:\t%g\n", prod*exp(large_num1));
  
  //fprintf(stdout,"(FB)/(Laprace):\t\t\t%f\n",fb/prod*exp(large_num2-large_num1));
  fprintf(stdout, "%10.6g ", fb/prod*exp(large_num2-large_num1) );
  print_vector(stdout, g, rank, "");

  return;
}

static void 
print_vector(FILE *fp, double *v, int length, char *str)
{
  int i;
  fprintf(fp, "%s\t", str);
  for (i=0; i<length; i++)
    fprintf(fp, "%10.6g ", v[i]);
  fprintf(fp, "\n");
}

static double 
gamma_function1(int k)
{
  if ( k == 1) 
    return 1.0;
  else 
    return (k-1)*gamma_function1(k-1);
}

static double 
gamma_function2(int k)
{
  if ( k == 1) 
    return sqrt(M_PI);
  else if ( k % 2 == 1)
    return (0.5*k-1)*gamma_function2(k-2);
  else 
    return 0.0;
}

static void 
comparison(void)
{
  double sigma_mu[6*dim];
  s  = sigma_mu;
  m  = sigma_mu +   dim;
  ss = sigma_mu + 2*dim;
  mm = sigma_mu + 3*dim;
  x  = sigma_mu + 4*dim;
  y  = sigma_mu + 5*dim;

  int i;
  for (i=0; i<dim; i++){
    s[i] = 1.0;
    m[i] = 0.0;
  }

  /* Set parameters */
  rank = 2 * dim;
  for ( i = 0; i < dim; i++){
    ss[i] = s[i] * s[i];
    mm[i] = m[i] * m[i];
    x[i] = -0.5 / (s[i]*s[i]);
    y[i] = m[i] / (s[i]*s[i]);    
  }
  /* initial value */
  double r = INITIAL_R;
  double f[rank];
  for (i = 0; i < dim; i++){
    f[i]     = y[i];
    f[i+dim] = 1.0;
  }

  /* set Log_C*/
  Log_C= (dim+1)*log(r);
  for (i = dim; i > 0; i -= 2) Log_C -= log(i);
  if (dim % 2) {
    Log_C += log(2.0) + 0.5*(dim-1)*log(2*M_PI);
  } else {
    Log_C += 0.5*dim*log(2*M_PI);
  }
  /* move r form 1e-6 to 1.0 */  
  runge_kutta_3(f, INITIAL_R, 1.0);
  r = 1.0;

  for ( i = 0 ; i < rank; i++)
    f[i] *= exp(-r*r*x[0]-r*y[0]);
  f[0] /= r;
  f[dim] /= r*r;

  /* move r form 1.0 to radius */
  runge_kutta_4(f, 1.0, radius);
  r = radius;

  /* Fisher--Bingham Integral */
  double fb;
  double large_num2=Log_C+r*r*x[0]+r*y[0];
  fb = r*r*f[dim];
  for ( i = 1; i < dim; i++)
    fb += f[dim+i];
  fb /= r*r;
  fb *= exp(large_num2);
  //fprintf(stdout, "Fisher-Bingham integral:\t%g\n", fb );

  double gamma;
  if (dim % 2 == 0)
    gamma = gamma_function1(dim/2);
  else 
    gamma = gamma_function2(dim);
  double fb2 = 2 * pow(M_PI, 0.5*dim) * pow(r,dim-1) * exp(-0.5*r*r) / gamma;
  fprintf(stdout, "%2d & %f & %4.2e\\\\\n", dim, fb, fb-fb2 );

  return;
}

void 
feller(void)
{
  int dim2 = dim;
  dim *= 2;
  //printf("dim=%d\n", dim);

  double sigma_mu[6*dim];
  s  = sigma_mu;
  m  = sigma_mu +   dim;
  ss = sigma_mu + 2*dim;
  mm = sigma_mu + 3*dim;
  x  = sigma_mu + 4*dim;
  y  = sigma_mu + 5*dim;

  int i;
  for (i=0; i<dim2; i++){
    s[2*i]   = s[2*i+1] = sqrt( 0.5/(i+1.0) );
  }
  for (i=0; i<dim; i++){
    m[i] = 0.0;
  }

  /* Set parameters */
  rank = 2 * dim + 1;
  for ( i = 0; i < dim; i++){
    ss[i] = s[i] * s[i];
    mm[i] = m[i] * m[i];
    x[i] = -0.5 / (s[i]*s[i]);
    y[i] = m[i] / (s[i]*s[i]);    
  }
  /*
    print_vector(stdout, s, dim, "  sigma:");
    print_vector(stdout, m, dim, "     mu:");
    print_vector(stdout, x, dim, "      x:");
    print_vector(stdout, y, dim, "      y:");
  */

  /* initial value */
  double f[rank];
  for (i = 0; i < dim; i++){
    f[i]     = y[i];
    f[i+dim] = 1.0;
  }
  f[rank-1] = 0.0;
  //  print_vector(stdout, f, rank, "      f:");
  double r = INITIAL_R;

  /* set Log_C*/
  Log_C= (dim+1)*log(r);
  for (i = 0; i < dim; i++)    Log_C -= log(s[i]);
  for (i = 0; i < dim; i++)    Log_C -= 0.5*mm[i]/ss[i];
  for (i = dim; i > 0; i -= 2) Log_C -= log(i);
  if (dim % 2) Log_C += 0.5 * (log(2.0) - log(M_PI));
  //  printf("Log_C: %.10g\n", Log_C);
  //  printf("    C: %.10g\n", exp(Log_C));

  if (radius < 1.0) {
    /* move r form 1e-6 to radius */  
    runge_kutta_1(f, INITIAL_R, radius);
    print_vector(stdout, f, rank, "f:\t");
    r = radius;
  } else {
    /* move r form 1e-6 to 1.0 */  
    runge_kutta_1(f, INITIAL_R, 1.0);
    r = 1.0;
    
    //  fprintf(stdout, "r=%f, C=%g\n", r, exp(Log_C));
    //  print_vector(stdout, f, rank, "f:\t");
    
    for ( i = 0 ; i < rank-1; i++)
      f[i] *= exp(-r*r*x[0]-r*y[0]);
    f[0] /= r;
    f[dim] /= r*r;
    //print_vector(stdout, f, rank, "f:\t");

    /* move r form 1.0 to radius */
    runge_kutta_2(f, 1.0, radius);
    r = radius;
  }

  //fprintf(stdout, "r=%f\n", r);
  //print_vector(stdout, f, rank, "f:\t");

  double fb1=f[rank-1];
  double fb2=pow(1.0-exp(-r*r), dim2);

  fprintf(stdout, "%2d & %f & %4.2e\\\\\n", dim, fb1, fb1-fb2 );
  return;
}

int
main(int argc, char *argv[])
{
  int target;
  sscanf(argv[1], "%d", &target);
  sscanf(argv[2], "%d", &dim);
  sscanf(argv[3], "%lf", &radius);

  double sigma_mu[6*dim];
  s  = sigma_mu;
  m  = sigma_mu +   dim;
  ss = sigma_mu + 2*dim;
  mm = sigma_mu + 3*dim;
  x  = sigma_mu + 4*dim;
  y  = sigma_mu + 5*dim;

  int i;
  char **p = argv+4;
  if ( target < 5){
  for ( i = 0; i < dim; i++)
    sscanf(*p++, "%lf", s+i);
  for ( i = 0; i < dim; i++)
    sscanf(*p++, "%lf", m+i);
  }

  switch (target) {
  case 1:
    ball_prob();
    break;
  case 2:
    ball_prob_graph();
    break;
  case 3:
    fisher_bingham();
    break;
  case 4:
    fisher_bingham_graph();
    break;
  case 5:
    hirotsu_1_1();
    break;
  case 6:
    hirotsu_1_2();
    break;
  case 7:
    hirotsu_2_1();
    break;
  case 8:
    hirotsu_2_2();
    break;
  case 9:
    anderson_darling();
    break;
  case 10:
    comparison();
    break;
  case 11:
    feller();
    break;
  default:
    break;
  }
  return 0;
}
