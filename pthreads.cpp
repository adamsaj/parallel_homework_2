#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <pthread.h>


namespace {
  pthread_barrier_t barr;
  
  FILE *fsave;
  int n;
  int n_threads;
  particle_t *particles;
  bin_t* bins;

  struct thread_data_t {
    int threadno;
    int threads;
  };
}

void thread_inner(int threadno, int threads);

void* thread_launch(void* data) {
  struct thread_data_t* typed = (struct thread_data_t*)data;
  thread_inner(typed->threadno, typed->threads);
  return NULL;
}

void thread_inner(int threadno, int threads) {
  const int block_size = NUM_BINS/threads + 1;
  const int lstart = threadno * block_size;
  const int lend = min(NUM_BINS, (threadno+1)*block_size);

  for( int step = 0; step < NSTEPS; step++ ) {
    // 
    // for each bin, locally compute forces on particles
    // need to look only at self and neighboring bins
    //
    for(int i = lstart; i < lend; i++) {
      int r = i/NUM_BINS_PER_SIDE;
      int c = i % NUM_BINS_PER_SIDE;
      compute_forces_for_box(bins, r, c); // O(#bins * c)
    }

    //need to do a barrier like thing here
    // Synchronization point
    int rc = pthread_barrier_wait(&barr);
    if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD) {
      printf("Could not wait on barrier\n");
      exit(-1);
    }

    if( threadno == 0 ) {
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

    rc = pthread_barrier_wait(&barr);
    if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD) {
      printf("Could not wait on barrier\n");
      exit(-1);
    }

  }
}

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
        printf( "-p <int> to set the number of threads\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    n = read_int( argc, argv, "-n", 1000 );
    n_threads = read_int( argc, argv, "-p", 2 );
    printf("threads: %d\n", n_threads);

    // Barrier initialization
    if(pthread_barrier_init(&barr, NULL, n_threads)) {
      printf("Could not create a barrier\n");
      return -1;
    }



    char *savename = read_string( argc, argv, "-o", NULL );
    
    fsave = savename ? fopen( savename, "w" ) : NULL;
    particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );

    bins = new bin_t[NUM_BINS];
    //    bin_t *bins = (bin_t*) malloc( NUM_BINS * sizeof(bin_t) );
    init_bins( bins ); // O(bins)

    // 
    // assign each particle to the appropriate bin
    // 
    for( int i = 0; i < n; i++ ) { // O(n)
      bin_particle(particles[i], bins);
    }
    
    pthread_t* thr = new pthread_t[n_threads];

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    //we pay the thread launch in the timed block, but not
    // worth the engineering overhead to fix
    for(int i = 1; i < n_threads; ++i) {
      thread_data_t* data = new thread_data_t();
      data->threadno = i;
      data->threads = n_threads;
      if(pthread_create(&(thr)[i], NULL, &thread_launch, (void*)data)) {
        printf("Could not create thread %d\n", i);
        return -1;
      }
    }

    thread_data_t* data = new thread_data_t();
    data->threadno = 0;
    data->threads = n_threads;
    thread_launch(data);

    for(int i = 1; i < n_threads; ++i) {
      if(pthread_join(thr[i], NULL)) {
        printf("Could not join thread %d\n", i);
        return -1;
      }
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
