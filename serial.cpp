#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );

    bin_t* bins = new bin_t[NUM_BINS];
    //    bin_t *bins = (bin_t*) malloc( NUM_BINS * sizeof(bin_t) );
    init_bins( bins ); // O(bins)

    // 
    // assign each particle to the appropriate bin
    // 
    for( int i = 0; i < n; i++ ) { // O(n)
      bin_particle(particles[i], bins);
    }
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
      // 
      // for each bin, locally compute forces on particles
      // need to look only at self and neighboring bins
      //

      //asm("int3;");
#if 1
      for(int r = 0; r < NUM_BINS_PER_SIDE; r++) {
        for(int c = 0; c < NUM_BINS_PER_SIDE; c++) {
          compute_forces_for_box(bins, r, c); // O(#bins * c)
        }
      }
#else
      //TODO: This is the horrible O(n^2) code, it's correct but slow 
      for( int i = 0; i < n; i++ ) {
        particles[i].ax = particles[i].ay = 0;
        for (int j = 0; j < n; j++ )
          apply_force( particles[i], particles[j] );
      }
#endif

      // for each particle, move O(n)
      for( int i = 0; i < n; i++ ) 
          move( particles[i] );

      // for each particle, rebin for next iteration
      init_bins( bins ); //clear points from last iteration O(bins)
      for( int i = 0; i < n; i++ ) { //rebin from scratch, O(n)
        bin_particle(particles[i], bins);
      }


      //
      //  save if necessary
      //
      if( fsave && (step%SAVEFREQ) == 0 )
        //if( fsave )
        save( fsave, n, particles );
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );

    // not really needed, but simply clear all pointers to particles
    init_bins( bins );   
    free( particles );
    delete[] ( bins );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
