#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <time.h>

#include "./utils/amut.h"
#include "./draw/amgr.h"
#include "./utils/ambs.h"

#define nb_valeurs 10000

int verbose = 2;

double f_star=-1000, f_max=-1000;
double brownian[nb_valeurs];

double max(double x, double y) {
	if (x>y) return x;
	return y;
}
double min(double x, double y) {
	if (x<y) return x;
	return y;
}
double absolute(double x) {
	if (x>0) return x;
	return -x;
}

/***************** Object function ****************/
double f_start=-5, f_end=5;
double f_center=0;
int devide = 10;
double eps = 1.0E-4;

double f(double value) {
    double x1,x2,x3;
    x1 = value*value - value*5.0 + 6;
    x2 = value*value + 1.0;
    x3 = x1/x2;
    return x3;
}

double f_lowerbound(double lb, double ub){
    double a1,a2,a3,a4,a5,a6,b1,b2;
    a1 = lb*lb - ub*5.0 + 6.0;
    a2 = ub*ub - lb*5.0 + 6.0;
    a3 = ub*ub - ub*5.0 + 6.0;
    a4 = -ub*5.0 + 6.0;
    a5 = max(lb*lb,ub*ub)+1.0;
    a6 = max(lb*lb,ub*ub) - lb*5.0 + 6.0;
    b1 = ub*ub + 1.0;
    b2 = lb*lb + 1.0;
    if (lb>=0 && a1>=0)
        return a1/b1;
    else if (lb>=0 && a1<0 && a2>0)
        return a1/b2;
    else if (lb>=0 && a2<=0)
        return a1/b2;
    else if (ub<=0 && a3>0)
        return a3/b2;
    else if (lb<0 && ub>0 && a4>=0)
        return a4/a5;
    else if (lb<0 && ub>0 && a4<0 && a6>0)
        return a4;
    else if (lb<0 && ub>0 && a6<=0)
        return a4;
}

double f_upperbound(double lb,double ub){
    double a1,a2,a3,a4,a5,a6,a7,b1,b2;
    a1 = lb*lb - ub*5.0 + 6.0;
    a2 = ub*ub - lb*5.0 + 6.0;
    a3 = ub*ub - ub*5.0 + 6.0;
    a4 = -ub*5.0 + 6.0;
    a5 = max(lb*lb,ub*ub)+1.0;
    a6 = max(lb*lb,ub*ub) - lb*5.0 + 6.0;
    a7 = lb*lb - lb*5.0 + 6.0;
    b1 = ub*ub + 1.0;
    b2 = lb*lb + 1.0;
    if (lb>=0 && a1>=0)
        return a2/b2;
    else if (lb>=0 && a1<0 && a2>0)
        return a2/b2;
    else if (lb>=0 && a2<=0)
        return a2/b1;
    else if (ub<=0 && a3>0)
        return a7/b1;
    else if (lb<0 && ub>0 && a4>=0)
        return a6;
    else if (lb<0 && ub>0 && a4<0 && a6>0)
        return a6;
    else if (lb<0 && ub>0 && a6<=0)
        return a6/a5;
}

/*double f(double x) {
    return  4*x*(1-x)*(0.75+0.25*(1-sqrt(absolute(sin(60.0*x)))));
}*/

/*double f(double value) {
   double x = brownian[(int) (value*nb_valeurs)];
   return x;
 }*/

/*double f(double x) {
   double y;
   if (x==0) return 1;
   y = 1-sqrt(x) + (-x*x +sqrt(x) )*(sin(1/pow(x,2))+1)/2;
   return y;
 }*/

/*double f(double x) {
  return max(3.6*x*(1-x),1-50*fabs(1-0.02-x));
}*/

/*double f(double x) {
  return x;
}*/

/*double f(double x) {
  int D=6, i;
  for (i=1; i<D+1; i++) {
	if (x < pow(2,-i))
		return 0.05 + ((double) (D-i))/ ((double) D);
	x = x-pow(2,-i);
  }
  return 1.05;
}*/

/*double f(double x) {
  int D=6, i;
  for (i=1; i<D+1; i++) {
	if (x < pow(3,-i))
		return 0.05 + ((double) (D-i))/ ((double) D);
	x = x-pow(3,-i);
  }
  return 1.05;
}*/

/***************** Profondeur maximale ******************/
unsigned int prof_max(unsigned int t) {
	return sqrt(t);
}

void draw_partition(dyv **x, dyv **x_min, dyv **x_max, ivec **leaf, ivec **new, dyv **y, int h_max)
{
  int i,h;
  for(h=0; h< h_max; h++) {
    for (i=0;i < dyv_size(x[h]); i++) {
      if (ivec_ref(leaf[h],i) == 1) {
	ag_set_pen_color(AG_BLACK);
	ag_line(dyv_ref(x_min[h],i),0.02,dyv_ref(x_min[h],i),-0.02);
	ag_line(dyv_ref(x_max[h],i),0.02,dyv_ref(x_max[h],i),-0.02);
	ag_set_pen_color(AG_BLUE);
	ag_dot(dyv_ref(x[h],i),dyv_ref(y[h],i));
      }
    }
  }
}

// void draw_function() {
// 	double x, y;
// 	for (x=0; x<1; x+= 1.0 / ((double) nb_valeurs)) {
// 		y = f(x);
// 		ag_pixel(x, y);
// 		if (y>f_star) f_star = y;
// 	}
// }

void draw_function(double x_min, double x_max) {
	double x, y;
	for (x=x_min; x<x_max; x+= (x_max-x_min) * 0.00002) {
		y = f(x);
		ag_pixel((x-x_min)/(x_max-x_min), y);
 		if (y>f_star) f_star = y;
	}
}

void draw_tree(dyv **x, dyv **x_min, dyv **x_max, ivec **leaf, ivec **new, dyv **y, int h_max)
{
  int i,h;
  ag_on("");
  set_ag_frame(0,0,1,h_max);
  ag_set_pen_color(AG_BLACK);

  for(h=0; h< h_max; h++) {
    for (i=0;i < dyv_size(x[h]); i++) {
      if (ivec_ref(leaf[h],i) == 1) {
	ag_set_pen_color(AG_BLACK);
	ag_line(dyv_ref(x[h],i),0,dyv_ref(x[h],i),h);
      }
    }
  }
}

void init_drawing()
{
        //ag_on("");
	ag_on("tree.eps");
  	set_ag_frame(-0.1, -0.1, 1.1, 1.1);
	ag_set_pen_color(AG_BLACK);
	ag_line(0, 0, 1, 0);
	ag_set_pen_color(AG_RED);
}


void save_partition(dyv **x, dyv **x_min, dyv **x_max, ivec **leaf, ivec **new, dyv **y, int h_max)
{
  ag_on("tree.eps");
  set_ag_frame(-0.1, -0.1, 1.1, 1.1);
  ag_set_pen_color(AG_BLACK);
  ag_line(0, 0, 1, 0);
  ag_set_pen_color(AG_RED);
  draw_function(0,1);
//  draw_partition(x, x_min, x_max, leaf, new, y, h_max);
  ag_off();
}

void algorithme(dyv **x, dyv **x_min, dyv **x_max, ivec **leaf, ivec **new, dyv **y, int h_max, int nb_iter)
{
  int i,h,n;
  int at_least_one = 1;
  for (n=0;(n<nb_iter) && (at_least_one==1);) {
    double y_max = -10000;

    if (verbose >1) {
      init_drawing();
      draw_function(0,1);
 // really_wait_for_key();
      draw_partition(x,x_min,x_max,leaf,new, y,h_max);
      wait_for_key();
    }

    at_least_one=0;
    for(h=0; (h< h_max-1) && (n<nb_iter); h++) {
      int i_max = -1;
      for (i=0;i < dyv_size(x[h]); i++) {
        if ((ivec_ref(leaf[h],i) == 1) && (ivec_ref(new[h],i)==0)) {
	  double y_hi = dyv_ref(y[h],i);
	  if (y_hi > y_max) {
	    y_max = y_hi;
	    i_max = i;
	  }
	}
      }
      if (i_max >= 0) { // On ouvre la feuille (h,i_max)
	  double x_g = (5 * dyv_ref(x_min[h],i_max) + dyv_ref(x_max[h],i_max))/6.0;
	  double x_d = (dyv_ref(x_min[h],i_max) + 5 * dyv_ref(x_max[h],i_max))/6.0;
	  double xx= dyv_ref(x[h],i_max), yy = dyv_ref(y[h],i_max);
	  //if (verbose >1) {
	  //  init_drawing();
	  //  draw_function();
	  //}
	  //draw_partition(x,x_min,x_max,leaf,new,y,h_max);
	  ag_set_pen_color(AG_RED);
	  ag_circle(xx,yy,5);
	  ag_set_pen_color(AG_BLACK);
	  ag_box(dyv_ref(x_min[h],i_max),-0.01, dyv_ref(x_max[h],i_max),0.01);
	  printf("Step %d: f(%g) = %g\n",n, xx, yy);
	  //wait_for_key();

	  if (yy>f_max) f_max = yy;

	  at_least_one = 1;
          ivec_set(leaf[h],i_max, 0); // La feuille devient un noeud
	  n++; // on incrémente le nombre d'itérations
	  // Noeud fils gauche:
	  add_to_dyv(x[h+1], x_g);
	  add_to_dyv(y[h+1], f(x_g));
	  add_to_dyv(x_min[h+1], dyv_ref(x_min[h],i_max));
	  add_to_dyv(x_max[h+1], (2*dyv_ref(x_min[h],i_max)+dyv_ref(x_max[h],i_max))/3.0);
	  add_to_ivec(leaf[h+1],1);
	  add_to_ivec(new[h+1],1);
	  // Noeud fils droit:
  	  add_to_dyv(x[h+1], x_d);
	  add_to_dyv(y[h+1], f(x_d));
	  add_to_dyv(x_min[h+1], (dyv_ref(x_min[h],i_max)+2*dyv_ref(x_max[h],i_max))/3.0);
	  add_to_dyv(x_max[h+1], dyv_ref(x_max[h],i_max));
	  add_to_ivec(leaf[h+1],1);
	  add_to_ivec(new[h+1],1);
	  // Noeud fils central:
  	  add_to_dyv(x[h+1], xx);
	  add_to_dyv(y[h+1], yy);
	  add_to_dyv(x_min[h+1], (2*dyv_ref(x_min[h],i_max)+dyv_ref(x_max[h],i_max))/3.0);
	  add_to_dyv(x_max[h+1], (dyv_ref(x_min[h],i_max)+2*dyv_ref(x_max[h],i_max))/3.0);
	  add_to_ivec(leaf[h+1],1);
	  add_to_ivec(new[h+1],1);
	  //printf("Feuille (%d,%d) x_min=%g, x=%g, x_max=%g, y=%g \n",h,i_max,dyv_ref(x_min[h],i_max),dyv_ref(x[h],i_max),dyv_ref(x_max[h],i_max),dyv_ref(y[h],i_max));
      }
    }
    for(h=0; h< h_max; h++) {
      for (i=0;i < dyv_size(x[h]); i++) {
        ivec_set(new[h],i,0);
      }
    }
    ag_off();
    wait_for_key();
  }
}


int main(int argc, char* argv[]) {
	int h, nb_iter, h_max, i;

	if(argc != 2) {
		printf("Usage: %s nbIterations\n",argv[0]);
		printf("Ex: soo 100 for 100 iterations\n");
		return 1;
	}
	printf("Hi man 0");
	nb_iter = (unsigned int) atoi(argv[1]);

	dyv *x, *x_min, *x_max;
	dyv *y, *lb, *ub;
	*x = mk_dyv(0);
	*x_min = mk_dyv(0);
	*x_max = mk_dyv(0);
	*y = mk_dyv(0);
	*lb = mk_dyv(0);
	*ub = mk_dyv(0);

    double m,l;
    m = (f_end-f_start)/devide;
    l = m/2;
	for(i=0;i<devide;++i){
	    add_to_dyv(*x_min, f_start+i*m);
	    add_to_dyv(*x, f_start+i*m+l);
        add_to_dyv(*x_max, f_start+i*m+m);
        add_to_dyv(*y, f(f_start+i*m+l));
        add_to_dyv(*lb, f_lowerbound(f_start+i*m,f_start+i*m+m));
        add_to_dyv(*lb, f_upperbound(f_start+i*m,f_start+i*m+m));
	}
}





/*
	{
	dyv *x[h_max], *x_min[h_max], *x_max[h_max];
	ivec *leaf[h_max], *new[h_max];
	dyv *y[h_max];

	for (h=0; h<h_max; h++) {
	  x[h] = mk_dyv(0);
	  x_min[h] = mk_dyv(0);
	  x_max[h] = mk_dyv(0);
	  leaf[h] = mk_ivec(0);
	  new[h] = mk_ivec(0);
	  y[h] = mk_dyv(0);
	}
	add_to_dyv(x_min[0], 0);
	add_to_dyv(x_max[0], 1);
	add_to_dyv(x[0], 0.5);
    add_to_ivec(leaf[0],1);
    add_to_ivec(new[0],0);
	//definition_brownian();
	add_to_dyv(y[0],f(0.5));

	if (verbose >1) {
	  init_drawing();
	  //draw_function(0,1);
	}
	algorithme(x,x_min,x_max,leaf,new,y,h_max,nb_iter);
	//printf("\nRegret: %g\n", f_star - f_max);
	printf("Max value: %.100f\n", f_max);
	printf("Regret: %.100f\n", 0.975599143811574975870826165191829204559326171875 - f_max);

	save_partition(x,x_min,x_max,leaf,new,y,h_max);

	wait_for_key();
	draw_tree(x,x_min,x_max,leaf,new,y,h_max);

	wait_for_key();

	// Libère toute les places mémoire
	for (h=0; h<h_max; h++) {
	  free_dyv(x[h]); free_dyv(x_min[h]); free_dyv(x_max[h]);
          free_ivec(leaf[h]); free_ivec(new[h]); free_dyv(y[h]);
	}
	}
	return 0;
}
*/
