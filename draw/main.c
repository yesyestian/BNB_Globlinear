#include "amut.h"
#include "amgr.h"
#include "ammarep.h"
#include "drac.h"

typedef struct gauss_info
{
  dyv *mu;
  dym *cov_inv;
  double cov_determinant;
} gauss_info;

double gauss_eval(dyv *x,dyv *mu,dym *cov_inv,double covdet)
{
  dyv *x_minus_mu = mk_dyv_subtract(x,mu);
  double quad_form = dym_xt_a_x_value(x_minus_mu,cov_inv);
  double multiplicand = 1 / sqrt(2 * PI * covdet);
  double result = multiplicand * exp(-quad_form / 2);
  free_dyv(x_minus_mu);
  return result;
}  
  
double gauss_height_fn(char *data,double x,double z)
{
  gauss_info *gi = (gauss_info *) data;
  dyv *xdyv = mk_dyv_2(x,z);
  double result = gauss_eval(xdyv,gi->mu,gi->cov_inv,gi->cov_determinant);
  free_dyv(xdyv);
  return result;
}

surgraph *mk_surgraph_from_2d_function(
    double (*height_function)(char *data, double x, double z),
    char *data,
    int num_contours,
    int grid_size,
    char *title,
    char *xlabel,
    char *zlabel,
    double x_low,
    double z_low,
    double x_high,
    double z_high
  );

void draw_2d_gaussian(dyv *mu,dym *cov,dyv *lo,dyv *hi)
{
  surgraph *sg;
  gauss_info gi[1];
  char buff[1000];

  gi -> mu = mu;
  gi -> cov_inv = mk_invert_dym(cov);
  gi -> cov_determinant = dym_determinant(cov);

  sprintf(buff,"mu = (%g,%g), cov=((%g,%g),(%g,%g))",
          dyv_ref(mu,0),dyv_ref(mu,1),
          dym_ref(cov,0,0),dym_ref(cov,0,1),
          dym_ref(cov,1,0),dym_ref(cov,1,1));

  sg = mk_surgraph_from_2d_function(gauss_height_fn,(char *)gi,
               30,30,buff,"x1","x2",dyv_ref(lo,0),dyv_ref(lo,1),
               dyv_ref(hi,0),dyv_ref(hi,1));

  ag_on("gauss.ps");
  render_surgraph(sg);
  ag_off();

  free_dym(gi->cov_inv);
  free_surgraph(sg);
}

void gauss_main(int argc,char *argv[])
{
  dyv *dmu = mk_dyv_2(1.0,2.0);
  dym *dcov = mk_dym_22(1.0,0.0,0.0,3.0);
  dyv *dlo = mk_dyv_2(-4.0,-4.0);
  dyv *dhi = mk_dyv_2(4.0,4.0);
  dyv *mu = mk_dyv_from_args("mu",argc,argv,dmu);
  dym *cov = mk_dym_from_args("cov",argc,argv,dcov);
  dyv *lo = mk_dyv_from_args("lo",argc,argv,dlo);
  dyv *hi = mk_dyv_from_args("hi",argc,argv,dhi);

  draw_2d_gaussian(mu,cov,lo,hi);

  free_dyv(dmu);
  free_dym(dcov);
  free_dyv(dlo);
  free_dyv(dhi);
  free_dyv(mu);
  free_dym(cov);
  free_dyv(lo);
  free_dyv(hi);

  wait_for_key();
  am_malloc_report();
}

void screen_main(int argc,char *argv[])
{
  int i;
  printf("Will do ag_window_shape(384,384);\n");
  wait_for_key();
  ag_window_shape(384,384);

  printf("Will do ag_on(test.ps);\n");
  wait_for_key();
  ag_on("test.ps");

  printf("Will do ag_line(10.0,10.0,256.0,502.0);\n");
  wait_for_key();
  ag_line(10.0,10.0,256.0,502.0);

  printf("Will do ag_off();\n");
  wait_for_key();
  ag_off();

  printf("Will do ag_on("");\n");
  wait_for_key();
  ag_on("");

  printf("Will do ag_line(10.0,10.0,502.0,256.0);\n");
  wait_for_key();
  ag_line(10.0,10.0,502.0,256.0);

  printf("Will do ag_screen_shape(600.0,600.0);\n");
  wait_for_key();
  ag_window_shape(600.0,600.0);

  printf("Will do ag_line(10.0,10.0,502.0,256.0);\n");
  wait_for_key();
  ag_line(10.0,10.0,502.0,256.0);

  printf("Will do ag_off();\n");
  wait_for_key();
  ag_off();

  printf("Will do ag_on(test2.ps);\n");
  wait_for_key();
  ag_on("test2.ps");

  printf("Will do ag_line(10.0,10.0,256.0,502.0);\n");
  wait_for_key();
  ag_line(10.0,10.0,256.0,502.0);

  for ( i = 2 ; i < 25 ; i++ )
  {
    if ( i == 12 ) ag_set_pen_color(AG_RED);
    ag_pixel((double)i,(double)(i*i));
  }

  printf("Will do ag_off();\n");
  wait_for_key();
  ag_off();

  printf("Done\n");
  wait_for_key();
}

void test_pats()
{
  ps_usepats();
  ag_on("pattest.ps");
  ag_set_pen_color(AG_BLACK);
  ag_disc(100,100,50);
  ag_set_pen_color(AG_WHITE);
  ag_disc(100,100,40);
  ag_set_pen_color(AG_RED);
  ag_disc(100,100,30);
  ag_set_pen_color(AG_WHITE);
  ag_disc(100,100,20);
  ag_set_pen_color(AG_BLACK);
  ag_disc(100,100,10);
  ag_rectangle(50,50,150,150);
  printf("This bullseye should show up in pattest.ps, with the red part\n");
  printf("turned into a pattern instead of a color\n");
  wait_for_key();
  ag_off();
}

int main(int argc,char *argv[])
{
  void screen_main(int argc,char *argv[]);
  test_pats();
  screen_main(argc,argv);
  return 0;
}
