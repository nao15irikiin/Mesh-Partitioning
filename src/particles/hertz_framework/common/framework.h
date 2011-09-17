#ifndef FRAMEWORK_H
#define FRAMEWORK_H

#ifdef GPU_TIMER
  #include "cuda_timer.h"
#else
  #include "simple_timer.h"
#endif

#include "unpickle.h"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

std::vector<SimpleTimer> one_time;
std::vector<SimpleTimer> per_iter;

/* return a random integer in [0..n) */
int rand_int(int n) {
  int limit = RAND_MAX - RAND_MAX % n;
  int rnd;

  do {
    rnd = random();
  } while (rnd >= limit);
  return rnd % n;

}

/* randomly shuffle the contact list (array of pairs) */
void shuffle_edges(int *edges, int nedge) {
  int i, j, n0, n1;

  for (i = nedge - 1; i > 0; i--) {
    j = rand_int(i + 1);
    n0 = edges[(j*2)];
    n1 = edges[(j*2)+1];
    edges[(j*2)]   = edges[(i*2)];
    edges[(j*2)+1] = edges[(i*2)+1];
    edges[(i*2)]   = n0;
    edges[(i*2)+1] = n1;
  }
}

/*
 * Run [num_iter] iterations of the hertz computation. [one_time] and [per_iter]
 * store timing results for one-time and per-iteration costs, respectively.
 */
extern void run(struct params *input, int num_iter);

void print_usage(std::string progname) {
  printf("Usage: %s <stepfile> [options]\n", progname.c_str());
  printf("Options:\n");
  printf("   -n arg     number of runs\n");
  printf("   -v         be verbose\n");
  printf("   -p arg     use partition file\n");
  printf("   -s arg     set seed for edge shuffle\n");
}

int main(int argc, char **argv) {

  // PARSE CMDLINE
  std::string progname(argv[0]);

  // mandatory arguments
  if (argc < 2) {
    print_usage(progname);
    return(1);
  }
  std::string step_filename(argv[1]);
  struct params *p = parse_file(step_filename);
  argc--;
  argv++;

  // optional arguments
  bool debug = false;
  bool verbose = false;
  int num_iter = 1000;
  std::string part_filename;
  long seed = -1;

  bool isShuffle=false;		// Added by Nao
  bool isSmartPart=false;	// Added by Nao
extern int OP_part_size;
  int c;
  while ((c = getopt (argc, argv, "hdvn:p:s:b:")) != -1) {
    switch (c) {
      case 'h':
        print_usage(progname);
        return 1;
      case 'd':
        debug = true;
        break;
      case 'v':
        verbose = true;
        break;
      case 'n':
        num_iter = atoi(optarg);
        break;
      case 'p':
        part_filename = optarg;
        parse_partition_file(p, part_filename);
        isSmartPart=true;	// Added by Nao
        break;
      case 's':
        seed = atol(optarg);
        srandom(seed);
        shuffle_edges(p->edge, p->nedge);
        isShuffle=true;		// Added by Nao
        break;
      case 'b':
printf("found -b in main!!!\n");
OP_part_size = atoi(optarg);
        break;
      case '?':
        if (optopt == 'n' || optopt == 'p' || optopt == 's')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
              "Unknown option character `\\x%x'.\n",
              optopt);
        return 1;
      default:
        abort ();
    }
  }

  // pass argument for expansion in run function 
  p->argc=argc;
  p->argv=argv;

  if (debug) {
    printf ("# Command-line parsing: step_filename=%s verbose=%d num_iter=%d part_filename=%s seed=%ld\n",
        step_filename.c_str(), verbose, num_iter, part_filename.c_str(), seed);
    for (int i=optind; i<argc; i++)
      printf ("# Non-option argument: %s\n", argv[i]);
  }

  if (verbose) {
    printf("# Program: %s\n", progname.c_str());
    printf("# Num Iterations: %d\n", num_iter);
    if (p->npartition > 1) {
      printf("# Partition: %s\n", part_filename.c_str());
      printf("# npartition: %d\n", p->npartition);
    }
    if (seed != -1) {
      printf("# Shuffle seed: %ld\n", seed);
    }
#ifdef GPU_TIMER
    printf("# GPU timer implementation\n");
#elif POSIX_TIMER
    printf("# POSIX timer implementation\n");
#else
    printf("# CPU timer implementation\n");
#endif
  }

  // RUN TEST
  run(p, num_iter);
  double one_time_total = 0.0f;
  double per_iter_total = 0.0f;
  for (int i=0; i<one_time.size(); i++) {
    one_time_total += one_time[i].total_time();
  }
  for (int i=0; i<per_iter.size(); i++) {
    per_iter_total += per_iter[i].total_time();
  }

  if (verbose) { //then print header
    printf("# nedge, total_one_time_cost (milliseconds), time_per_iteration");
    for (int i=0; i<one_time.size(); i++) {
      printf(", [%s]", one_time[i].get_name().c_str());
    }
    for (int i=0; i<per_iter.size(); i++) {
      printf(", %s", per_iter[i].get_name().c_str());
    }
    if (seed != -1) {
      printf(", seed");
    }
    printf("\n");
  }

  // print runtime data
  printf("%d, %f, %f", p->nedge, one_time_total, per_iter_total / (double) num_iter);
  for (int i=0; i<one_time.size(); i++) {
    printf(", %f", one_time[i].total_time());
  }
  for (int i=0; i<per_iter.size(); i++) {
    printf(", %f", per_iter[i].total_time() / (double) num_iter);
  }
  if (seed != -1) {
    printf(", %ld", seed);
  }
  printf("\n");

  //For logging:Nao
  char *per_path;
  per_path = getenv("HERTZ_RSLT_PATH");

  char perfile[256];
  if(per_path == NULL)
    strcpy(perfile, "./hertz-output");
  else
    strcpy(perfile, per_path);

  strcat(perfile, "/perresult.dat");

  bool isFileExist=true;
  if(access(perfile, R_OK)==-1) isFileExist=false;

  FILE *fp_per;
  if ( (fp_per = fopen(perfile,"a")) == NULL) {
    printf("can't generate/find file %s\n", perfile); exit(-1);
  }

  if(isFileExist==false) {
    fprintf(fp_per, "<step file>, <is shuffle>, <partition type>, <part file>, <num of blocks>, <block size>, <num of colors>, <num of threadcolors>, <max block size>, <shared mem>");
    fprintf(fp_per, ", edges, total_one_time_cost (milliseconds), time_per_iteration");
    for (int i=0; i<one_time.size(); i++) {
      fprintf(fp_per, ", [%s]", one_time[i].get_name().c_str());
    }
    for (int i=0; i<per_iter.size(); i++) {
      fprintf(fp_per, ", %s", per_iter[i].get_name().c_str());
    }
    fprintf(fp_per, ", <time of exec>\n");
  }

  time_t timer;
  time(&timer);
  tm *ct = localtime(&timer);

  op_plan* myplan=op_plan_getByName("res");
  if(myplan == NULL) {printf("Can not find specified plan"); exit(-1);}

  // Print "res" data
  fprintf(fp_per, "%s, %s, %s, %s, %d, %d, %d, %d, %d, %.2fKB",
			step_filename.c_str(),
          (isShuffle)?"true":"false",
          (isSmartPart)?"smart_partition":"sequential_partition",
          (part_filename.empty())?"None":part_filename.c_str(),
          myplan->nblocks,
          BSIZE,
          myplan->ncolors,
          myplan->max_nthrcol,
          myplan->maxbsize,
          myplan->nshared/1024.0f);

  fprintf(fp_per, ", %d, %f, %f", p->nedge, one_time_total, per_iter_total / (double) num_iter);
  for (int i=0; i<one_time.size(); i++) {
    fprintf(fp_per, ", %f", one_time[i].total_time());
  }
  for (int i=0; i<per_iter.size(); i++) {
    fprintf(fp_per, ", %f", per_iter[i].total_time() / (double) num_iter);
  }

  fprintf(fp_per, ", %s",asctime(ct));

  fclose(fp_per);

  delete_params(p);
  return 0;

}
#endif
