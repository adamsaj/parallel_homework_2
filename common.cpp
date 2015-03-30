#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"


/// The size of one side of the square space
double size;

/// The size of each given bin.
double bin_size;
/// integer number of bins (with some slop)
int bins_per_side;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size( int n )
{
    size = sqrt( density * n );




    bin_size = sqrt(density*5);
    bins_per_side = floor(size/bin_size)+1;
    int num_bins = bins_per_side*bins_per_side;

    assert( bin_size > 2*cutoff );

    printf("%d, %d, %g\n", num_bins, bins_per_side, bin_size);
    printf("%g\n", size);
    printf("%d\n", n);
}

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p )
{
  //srand48( time( NULL ) );
    srand48( 1000 );    
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;
    
    for( int i = 0; i < n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];
        
        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);

        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
        p[i].globalID=i;

    }
    free( shuffle );
}

//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor )
{

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff )
        return;
    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );

    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }
}

//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    fprintf(f, "----\n");
    for( int i = 0; i < n; i++ )
        fprintf( f, "%10.8f %10.8f\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}

bin_t::bin_t() 
  : num_particles(0)
{}

//
// Added for homework
// 

void init_bins(bin_t* bins) {
  for( int i = 0; i < NUM_BINS; i++) {
    bins[i].num_particles = 0;
    bins[i].particles.clear();
  }
}

void bin_particle(particle_t &particle, bin_t* bins) {
  int col = floor(particle.x/BIN_SIZE);
  int row = floor(particle.y/BIN_SIZE);

  int index = row*NUM_BINS_PER_SIDE + col;
  bins[index].particles.push_back( &particle );
  bins[index].num_particles++;
}

///Compute all of the forces for the particles local to a specific box
/// i is the box index
void compute_forces_for_box(bin_t* boxes, int r, int c) {
  //Computing forces between all points in the box and all points in the neighbouring boxes:
  int i = r * NUM_BINS_PER_SIDE + c;
  int box_size = NUM_BINS_PER_SIDE;

  for(int p = 0; p < boxes[i].num_particles; p++) {

    (*boxes[i].particles[p]).ax = 0;
    (*boxes[i].particles[p]).ay = 0;

#define COMPARE_TO_BIN(index) \
    for(int p2 = 0; p2 < boxes[index].num_particles; p2++) \
      apply_force( *boxes[i].particles[p], *boxes[index].particles[p2] );

    //all particles in this bin 
    COMPARE_TO_BIN(i);

    //For the box to the left:
    if( c-1 >= 0 )
      COMPARE_TO_BIN(i-1);
    //For the box to the right:
    if( c+1 < NUM_BINS_PER_SIDE )
      COMPARE_TO_BIN(i+1);
    //For the box to the top:
    if (r-1 >= 0)
      COMPARE_TO_BIN(i-box_size);
    //For the box to the bottom:
    if (r+1 < NUM_BINS_PER_SIDE)
      COMPARE_TO_BIN(i+box_size);

      //For the box int the upper left-hand corner:
    if ((c - 1 >= 0) && (r-1 >= 0))
      COMPARE_TO_BIN(i -box_size -1 );
      //For the box int the upper right-hand corner:
    if ((c + 1 < NUM_BINS_PER_SIDE) && (r-1 >= 0))
      COMPARE_TO_BIN(i-box_size+1);
      //For the box int the lower left-hand corner:
    if ((c - 1 >= 0) && (r+1 < NUM_BINS_PER_SIDE))
      COMPARE_TO_BIN(i+box_size-1);
    //For the box int the lower right-hand corner:
    if ((c + 1 < NUM_BINS_PER_SIDE) && (r+1 < NUM_BINS_PER_SIDE))
      COMPARE_TO_BIN(i+box_size+1);
#undef COMPARE_TO_BIN
  }
}




mpi_bin_t::mpi_bin_t() 
  : num_particles(0)
{}


void mpi_init_bins(mpi_bin_t* bins) {
  for( int i = 0; i < NUM_BINS; i++) {
    bins[i].num_particles = 0;
    bins[i].particles.clear();
  }
}



void mpi_bin_particle(particle_t &particle, mpi_bin_t* bins) {
  int col = floor(particle.x/BIN_SIZE);
  int row = floor(particle.y/BIN_SIZE);

  int index = row*NUM_BINS_PER_SIDE + col;
  bins[index].particles.push_back( particle );
  bins[index].num_particles++;
}

///Compute all of the forces for the particles local to a specific box
/// i is the box index
void mpi_compute_forces_for_box(mpi_bin_t* boxes, int r, int c) {
  //Computing forces between all points in the box and all points in the neighbouring boxes:
  int i = r * NUM_BINS_PER_SIDE + c;
  int box_size = NUM_BINS_PER_SIDE;

  for(int p = 0; p < boxes[i].particles.size(); p++) {

    (boxes[i].particles[p]).ax = 0;
    (boxes[i].particles[p]).ay = 0;

#define COMPARE_TO_BIN(index) \
    for(int p2 = 0; p2 < boxes[index].particles.size(); p2++)           \
      apply_force( boxes[i].particles[p], boxes[index].particles[p2] );

    //all particles in this bin 
    COMPARE_TO_BIN(i);

    //For the box to the left:
    if( c-1 >= 0 )
      COMPARE_TO_BIN(i-1);
    //For the box to the right:
    if( c+1 < NUM_BINS_PER_SIDE )
      COMPARE_TO_BIN(i+1);
    //For the box to the top:
    if (r-1 >= 0)
      COMPARE_TO_BIN(i-box_size);
    //For the box to the bottom:
    if (r+1 < NUM_BINS_PER_SIDE)
      COMPARE_TO_BIN(i+box_size);

      //For the box int the upper left-hand corner:
    if ((c - 1 >= 0) && (r-1 >= 0))
      COMPARE_TO_BIN(i -box_size -1 );
      //For the box int the upper right-hand corner:
    if ((c + 1 < NUM_BINS_PER_SIDE) && (r-1 >= 0))
      COMPARE_TO_BIN(i-box_size+1);
      //For the box int the lower left-hand corner:
    if ((c - 1 >= 0) && (r+1 < NUM_BINS_PER_SIDE))
      COMPARE_TO_BIN(i+box_size-1);
    //For the box int the lower right-hand corner:
    if ((c + 1 < NUM_BINS_PER_SIDE) && (r+1 < NUM_BINS_PER_SIDE))
      COMPARE_TO_BIN(i+box_size+1);
#undef COMPARE_TO_BIN
  }
}

void mpi_get_valid_neighbor_blocks(int r, int c, std::vector<int>& out) {
  int i = r * NUM_BINS_PER_SIDE + c;
  int box_size = NUM_BINS_PER_SIDE;

  out.clear();
  //For the box to the left:
  if( c-1 >= 0 )
    out.push_back(i-1);
  //For the box to the right:
  if( c+1 < NUM_BINS_PER_SIDE )
    out.push_back(i+1);
  //For the box to the top:
  if (r-1 >= 0)
    out.push_back(i-box_size);
  //For the box to the bottom:
  if (r+1 < NUM_BINS_PER_SIDE)
    out.push_back(i+box_size);
  
  //For the box int the upper left-hand corner:
  if ((c - 1 >= 0) && (r-1 >= 0))
    out.push_back(i -box_size -1 );
  //For the box int the upper right-hand corner:
  if ((c + 1 < NUM_BINS_PER_SIDE) && (r-1 >= 0))
    out.push_back(i-box_size+1);
  //For the box int the lower left-hand corner:
  if ((c - 1 >= 0) && (r+1 < NUM_BINS_PER_SIDE))
    out.push_back(i+box_size-1);
  //For the box int the lower right-hand corner:
  if ((c + 1 < NUM_BINS_PER_SIDE) && (r+1 < NUM_BINS_PER_SIDE))
    out.push_back(i+box_size+1);
}            

