#include "amut.h"
#include "hash.h"

void mdp()
{
  double r1 = 4.0;
  double r2 = -1.0;
  double r4 = 2.0;
  double v1 = 0.0;
  double v2 = 0.0;
  double v3 = 0.0;
  double v4 = 0.0;
  while ( TRUE )
  {
    double q1u = r1 + v2;
    double q1d = r1 + 0.5 * v4;
    double q2u = r2 + 0.3333333333333 * (v1 + v2);
    double q2d = r2 + 0.5 * v4;
    double q4u = r4 + v2;
    double q4d = r4 + 0.5 * (v4 + v2);

    printf("%g %g %g %g\n",v1,v2,v3,v4);

    v1 = real_max(q1u,q1d);
    v2 = real_max(q2u,q2d);
    v4 = real_max(q4u,q4d);

    wait_for_key();
  }
}

dyv *mk_bootstrap_dyv(dyv *x)
{
  int size = dyv_size(x);
  dyv *b = mk_dyv(size);
  int i;
  for ( i = 0 ; i < size ; i++ )
     dyv_set(b,i,dyv_ref(x,int_random(size)));
  return b;
}

dyv *mk_random_dyv_in_unit_cube(int size)
{
  dyv *x = mk_dyv(size);
  int i;
  for ( i = 0 ; i < size ; i++ )
    dyv_set(x,i,range_random(0.0,1.0));
  return x;
}

int num_bootstraps_above_zero(int num_cats,int num_bootstraps)
{
  dyv *x1 = mk_random_dyv_in_unit_cube(num_cats);
  dyv *x2 = mk_random_dyv_in_unit_cube(num_cats);
  dyv *diff = mk_dyv_subtract(x1,x2);
  int i;
  int result = 0;

  for ( i = 0 ; i < num_bootstraps ; i++ )
  {
    dyv *boot = mk_bootstrap_dyv(diff);
    if ( dyv_mean(boot) > 0.1 )
      result += 1;
    free_dyv(boot);
  }

  free_dyv(diff);
  free_dyv(x1);
  free_dyv(x2);
  return result;
}

void boot_main(int argc,char *argv[])
{
  int num_trials = int_from_args("num_trials",argc,argv,300);
  int num_bootstraps = int_from_args("num_bootstraps",argc,argv,5000);
  int num_cats = int_from_args("num_cats",argc,argv,200);
  int i;
  for ( i = 0 ; i < num_trials ; i++ )
    printf("%9.6f\n",num_bootstraps_above_zero(num_cats,num_bootstraps) /
                     (double) num_bootstraps);
}

int main(int argc, char* argv[])
{
  printf("sqrt(2 * PI) = %20.18f\n",sqrt(2 * PI));

  if ((argc > 1) && (eq_string(argv[1],"boot")))
    boot_main(argc,argv);
  else if ((argc > 1) && (eq_string(argv[1],"mdp")))
    mdp();
  else if ((argc > 1) && (eq_string(argv[1],"hash_test")))
  {
    char s[200];
    hash_table *ht = mk_empty_hash_table(vint_copy,vint_free);
    string_array *sa;
    int done = 0;
    while (!done)
    {
      printf("> ");
      fflush(stdout);
      gets(s);
      sa = mk_broken_string(s);
      if (string_array_size(sa) < 1) continue;
      if (caseless_eq_string(string_array_ref(sa,0),"quit")) 
      {
	done = 1;
      }
      else if (caseless_eq_string(string_array_ref(sa,0),"lookup") &&
	       (string_array_size(sa) == 2))
      {
	int *a = (int *)lookup_hash_entry(ht,string_array_ref(sa,1));
	if (!a) printf("Not found in table.\n");
	else printf("Found %d.\n",*a);
      }
      else if (caseless_eq_string(string_array_ref(sa,0),"add") &&
	       (string_array_size(sa) == 3))
      {
	int a = atoi(string_array_ref(sa,2));
	add_hash_entry(ht,string_array_ref(sa,1),(void *)(&a));
      }
      else if (caseless_eq_string(string_array_ref(sa,0),"delete") &&
	       (string_array_size(sa) == 2))
      {
	delete_hash_entry(ht,string_array_ref(sa,1));
      }
      else
      {
	printf("Choices are:\n");
	printf("  quit\n");
	printf("  lookup key\n");
	printf("  add key number\n");
	printf("  delete key\n");
      }
      free_string_array(sa);
    }
    free_hash_table(ht);
    am_malloc_report();
  }
  if ((argc > 1) && (eq_string(argv[1],"string_sort_test")))
  {
    char s[200];
    string_array *sa;
    int done = 0;
    while (!done)
    {
      printf("> ");
      fflush(stdout);
      gets(s);
      sa = mk_broken_string(s);
      if (string_array_size(sa) < 1) continue;
      if (caseless_eq_string(string_array_ref(sa,0),"quit")) 
      {
	done = 1;
      }
      else
      {
	sort_string_array(sa);
	fprintf_string_array_contents_on_one_line(stdout,sa);
	fprintf(stdout,"\n");
      }
    }
  }
  else
  {
    printf("Nothing to see here, move along, move along...\n");
  }
  return 0;
}
