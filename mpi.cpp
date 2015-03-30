#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include <algorithm>
#include <vector>
#include <map>

//#define tracef printf
#define tracef //

#define ifdebug if(true) 
//#define ifdebug if(false) 


struct particle_key_less {
    bool operator ()(particle_t const& a, particle_t const& b) const {
      return a.globalID < b.globalID;
    }
};

/// Only meaningful when running on a single processor for debugging purposes
void assert_global_particle_count(mpi_bin_t* bins, int n, int rank, int n_proc, int my_bins_start, int my_bins_end) {
#if 1
  std::vector<particle_t> local;
  for(int r = 0; r < NUM_BINS_PER_SIDE; r++) {
    for(int c = 0; c < NUM_BINS_PER_SIDE; c++) {
      int i = r * NUM_BINS_PER_SIDE + c;
      
      if( my_bins_start <= i && i < my_bins_end ) {
        local.insert( local.end(), bins[i].particles.begin(), bins[i].particles.end() );
      }
    }
  }
  
  int lsize = local.size();
  assert( lsize <= n );
  //printf("lsize: %d\n", lsize);
  
  //get the sizes on the root
  int* sizes = new int[n_proc];   
  for(int i = 0; i < n_proc; i++) {
    sizes[i] = i; //reg bad data
  }
  MPI_Gather(&lsize, 1, MPI_INT, 
             sizes, 1, MPI_INT, 
             0, MPI_COMM_WORLD);

  //master proc will check and abort if needed
  if( rank == 0 ) {
    int gsize = 0;
    for(int i = 0; i < n_proc; i++) {
      gsize += sizes[i];
      //printf("proc %d has %d elements\n", i, sizes[i]);
    }
    if( n != gsize ) {
      printf("%d ==? %d\n", n, gsize);
      printf("\n\n\n");
      assert (gsize == n);
    }
  }
#else
  int result = 0;
  for( int i = 0; i < NUM_BINS; i++) {
    result += bins[i].particles.size();
  }
  assert (result == n);
#endif
}


//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    //
    //  process command line parameters
    //
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
    
    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    particle_t *particles = new particle_t[n];
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 7, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );


    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    set_size( n );
    if( rank == 0 )
        init_particles( n, particles );

    // Give all processors copies of the particles array (for the moment)
    MPI_Bcast(particles, n, PARTICLE, 0, MPI_COMM_WORLD);


    // Each processor fills in its own bin and a copy of ajacent ones
    // everything else is uninitialized and is kept around to make the code
    // copying easier from the serial version.
    mpi_bin_t* bins = new mpi_bin_t[NUM_BINS];
    mpi_init_bins( bins ); // O(bins) // TODO

    // 
    // assign each particle to the appropriate bin
    // done for all procs, all particles, we'll only
    // keep some of them up to date though
    // 
    for( int i = 0; i < n; i++ ) { // O(n)
      mpi_bin_particle(particles[i], bins); //TODO
    }

    delete[] particles;
    particles = NULL;

    // Figure out who owns what bin
    int bins_per_proc = NUM_BINS/n_proc + 1;
    int my_bins_start = bins_per_proc * rank;
    int my_bins_end = min(NUM_BINS, bins_per_proc * (rank +1) );
    //ex: rank 0 owns [0, bins_per_proc]

    // NOTE: Storing all bins on all procs is a bit of a hack.  We don't
    // update them, but most of the storage is wasted.


        //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ ) {
      
      //
      //  compute local forces 
      //
      for(int r = 0; r < NUM_BINS_PER_SIDE; r++) {
        for(int c = 0; c < NUM_BINS_PER_SIDE; c++) {
          int i = r * NUM_BINS_PER_SIDE + c;
          if( my_bins_start <= i && i < my_bins_end )
            mpi_compute_forces_for_box(bins, r, c); // O(#bins * c)
        }
      }
        
      //
      //  move local particles (by updating their properties, not rebinning)
      //
      for(int r = 0; r < NUM_BINS_PER_SIDE; r++) {
        for(int c = 0; c < NUM_BINS_PER_SIDE; c++) {
          int i = r * NUM_BINS_PER_SIDE + c;
          if( my_bins_start <= i && i < my_bins_end )
            for(int j = 0; j < bins[i].particles.size(); j++) {
              move( bins[i].particles[j] );
            }
        }
      }

#if 1
      // Now we need to exchange any updated points with our neighbors
      // We actually have to do this twice with possibly rebinning in
      // between - this catches the case where our neighbor has something
      // moved into them by their neighbor on the far side.
      for(int excnt = 0; excnt < 2; excnt++) {
        std::vector<MPI_Request> requests;

        //keep a copy of our old data (being sent) indexed by
        // bin number until all of the ongoing sends have completed.
        std::map<int, std::vector<particle_t> > scratch;

        for(int r = 0; r < NUM_BINS_PER_SIDE; r++) {
          for(int c = 0; c < NUM_BINS_PER_SIDE; c++) {
            int i = r * NUM_BINS_PER_SIDE + c;

            if( my_bins_start <= i && i < my_bins_end ) {
              //figure out the bin numbers of our neighbors (valid ones only)
              std::vector<int> neighbors;
              mpi_get_valid_neighbor_blocks(r, c, neighbors);

              //capture a copy until after all sends complete
              scratch[i] = bins[i].particles;

              //Send our content to the neighbors
              for(int j = 0; j < neighbors.size(); j++) {
                int owner = neighbors[j]/bins_per_proc;
         
                requests.push_back( MPI_Request() );
                
                //do a non-blocking send from that scratch space
                
                //trick: use bin number as tag 
                MPI_Isend(&(scratch[i][0]), scratch[i].size(), 
                          PARTICLE, owner, i, 
                          MPI_COMM_WORLD, &(requests.back()));

                tracef("Proc %d, sending bin %d with %d particles to proc %d\n",
                       rank, i, (int)scratch[i].size(), owner); 
              }
            }
          }
        }

        //For debugging only
        ifdebug
          assert_global_particle_count(bins, n, 
                                       rank, n_proc, 
                                       my_bins_start, my_bins_end);

        // MUST wait for all sends to complete before writing back into their
        // buffers again.
        std::vector<MPI_Status> junk;
        junk.resize(requests.size());
        //MPI_Waitall(requests.size(), &(requests[0]), &(junk[0]));

        //about to block on receives, must do all sends first
        // otherwise, we might write a buffer before reading it
        MPI_Barrier(MPI_COMM_WORLD);
        tracef("MPI_BARRIER\n");
        
        //do all the recieves in bulk (could make these irecv and remove
        // the barrier (todo)
        for(int r = 0; r < NUM_BINS_PER_SIDE; r++) {
          for(int c = 0; c < NUM_BINS_PER_SIDE; c++) {
            int i = r * NUM_BINS_PER_SIDE + c;

            if( my_bins_start <= i && i < my_bins_end ) {
              //figure out the bin numbers of our neighbors (valid ones only)
              std::vector<int> neighbors;
              mpi_get_valid_neighbor_blocks(r, c, neighbors);

              // wait for each neighbor to send us content
              for(int j = 0; j < neighbors.size(); j++) {
                int owner = neighbors[j]/bins_per_proc;

                //block waiting for the message (to figure out how large it is)
                MPI_Status status;
                MPI_Probe(owner, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                int len = 0;
                MPI_Get_count(&status, PARTICLE, &len);
                assert( len >= 0 );
                assert( len <= n );

                //trick: use bin number as tag
                int bin = status.MPI_TAG;
                assert( bin < NUM_BINS );
                assert( 0 <= bin );
                
                bins[bin].particles.resize(len);
                bins[bin].num_particles = len;

                //assert bin in neighbors

                MPI_Recv(&(bins[bin].particles[0]), bins[bin].particles.size(), 
                         PARTICLE, owner, bin, 
                         MPI_COMM_WORLD, &status);
                tracef("Proc %d, receiving bin %d with %d particles from proc %d\n",
                       rank, bin, len, owner); 
              }
            }
          }
        }
        //For debugging only
        ifdebug
          assert_global_particle_count(bins, n, 
                                       rank, n_proc, 
                                       my_bins_start, my_bins_end);

        // Since we just swapped info with all our neighbors, we need to rebin
        // all of our particles.
        std::vector<int> local_bins;
        for(int r = 0; r < NUM_BINS_PER_SIDE; r++) {
          for(int c = 0; c < NUM_BINS_PER_SIDE; c++) {
            int i = r * NUM_BINS_PER_SIDE + c;
              
            if( my_bins_start <= i && i < my_bins_end ) {
              std::vector<int> neighbors;
              mpi_get_valid_neighbor_blocks(r, c, neighbors);
              local_bins.insert( local_bins.end(), neighbors.begin(), neighbors.end() );
            }
          }
        }
        //remove duplicates to avoid double counting
        std::sort(local_bins.begin(), local_bins.end() );
        local_bins.erase(std::unique(local_bins.begin(), 
                                     local_bins.end()), 
                         local_bins.end());
        
        
        //grab all our locally tracked particles so we can rebin them
        std::vector<particle_t> local;
        for(int bi = 0; bi < local_bins.size(); bi++) {
          int i = local_bins[bi];
          local.insert( local.end(), 
                        bins[i].particles.begin(), 
                        bins[i].particles.end() );
          bins[i].particles.clear();
          bins[i].num_particles = 0;
        }
          
        for(int i = 0; i < local.size(); i++) {
          //if it's moved out of our area, essentially ignore it
          //otherwise replace it with our new info
          mpi_bin_particle(local[i], bins); 
        }
        //everything local should be correctly binned now
#if 1
        //debugging only - normally, this could overlap multiple iterations 
        MPI_Barrier(MPI_COMM_WORLD);
        int flag;
        //        std::vector<MPI_Status> junk;
        junk.resize(requests.size());
        MPI_Testall(requests.size(), &(requests[0]), &flag, &(junk[0]));
        assert( flag ); //otherwise recv missed. Ouch.
#endif
      }
#endif
      //For debugging only
      ifdebug
        assert_global_particle_count(bins, n, 
                                     rank, n_proc, 
                                     my_bins_start, my_bins_end);
      

      //
      //  save current step if necessary
      //
      if( savename && step%SAVEFREQ==0 ) {
        //combine all the local bins per proc
        std::vector<particle_t> local;
        for(int r = 0; r < NUM_BINS_PER_SIDE; r++) {
          for(int c = 0; c < NUM_BINS_PER_SIDE; c++) {
            int i = r * NUM_BINS_PER_SIDE + c;
          
            if( my_bins_start <= i && i < my_bins_end ) {
              local.insert( local.end(), bins[i].particles.begin(), bins[i].particles.end() );
              tracef("local bin %d of size %d\n", i, bins[i].particles.size());
            }
          }
        }

        int lsize = local.size();
        assert( lsize <= n );
      
        //get the sizes on the root
        int* sizes = new int[n_proc];   
        for(int i = 0; i < n_proc; i++) {
          sizes[i] = i; //reg bad data
        }
        MPI_Gather(&lsize, 1, MPI_INT, 
                   sizes, 1, MPI_INT, 
                   0, MPI_COMM_WORLD);
      
        //read back the data into a global buffer
      
        int gsize = 0;
        if( rank == 0 ) {
          for(int i = 0; i < n_proc; i++) {
            gsize += sizes[i];
            printf("proc %d, size %d\n", i, sizes[i]);
          }
          printf("%d == %d?\n", n, gsize);
          assert( n == gsize );
        }
        int* displ = NULL;
        if(rank == 0) {
          displ = new int[n_proc];
          displ[0] = 0;
          for(int i = 1; i < n_proc; i++) {
            displ[i] = displ[i-1]+sizes[i-1];
          }
        }
        std::vector<particle_t> global;
        if( rank == 0 )
          global.resize(gsize);

        MPI_Gatherv(&(local[0]), lsize, PARTICLE, 
                    &(global[0]), sizes, displ,
                    PARTICLE, 0, MPI_COMM_WORLD);
        local.clear();

        if( rank == 0 ) {
        
          //place them in order
          std::sort(global.begin(), global.end(), particle_key_less());
          
          for(int i = 1; i < global.size(); i++) {
            assert( global[i-1].globalID < global[i].globalID );
          }

          assert( fsave );

          save( fsave, n, &global[0] );
        }        
      }

      MPI_Barrier(MPI_COMM_WORLD);
    }
    simulation_time = read_timer( ) - simulation_time;
    
    if( rank == 0 )
        printf( "n = %d, n_procs = %d, simulation time = %g s\n", n, n_proc, simulation_time );
    
    //
    //  release resources
    //
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
