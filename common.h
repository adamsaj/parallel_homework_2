#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

#include <vector>

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 100;
const int SAVEFREQ = 10;

//
// particle data structure
//
typedef struct 
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
  int globalID;
} particle_t;

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
void init_particles( int n, particle_t *p );
void apply_force( particle_t &particle, particle_t &neighbor );
void move( particle_t &p );

//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

//
//  added for homework assignment
// 
//
// particle data structure
//
struct bin_t
{
  // This array holds the particles in a given bin
  // For simplicity, we only account for 10 particles
  // in a given bin.  Otherwise, we fall back to the slow code.
  // This is so we don't have to reallocate each time.
  std::vector<particle_t*> particles;
  int num_particles; //if 11 > then overflow and fall back

  bin_t();
};

/// The size of each given bin.
extern double bin_size;
/// integer number of bins (with some slop)
extern int bins_per_side;


#define BIN_SIZE bin_size
#define NUM_BINS_PER_SIDE bins_per_side
#define NUM_BINS (NUM_BINS_PER_SIDE*NUM_BINS_PER_SIDE)

/// Clear all the bins and remove any points from each
void init_bins(bin_t* bins);

/// Assign the given particle to the appropriate bin
void bin_particle(particle_t &particle, bin_t* bins);

void compute_forces_for_box(bin_t* boxes, int r, int c);


struct mpi_bin_t
{
  // Each particle lives in it's own local array only,
  // this makes movement more complicated, but removes the
  // need to update the whole particle array on every iteration
  std::vector<particle_t> particles;
  int num_particles; //for debugging purposes only

  mpi_bin_t();
};

/// Clear all the bins and remove any points from each
void mpi_init_bins(mpi_bin_t* bins);

/// Assign the given particle to the appropriate bin
void mpi_bin_particle(particle_t &particle, mpi_bin_t* bins);

void mpi_compute_forces_for_box(mpi_bin_t* boxes, int r, int c);

void mpi_get_valid_neighbor_blocks(int r, int c, std::vector<int>& out);
#endif
