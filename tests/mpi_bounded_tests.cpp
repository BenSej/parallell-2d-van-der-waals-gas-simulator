#include <random>

#include <mpi.h>

#include "gassim2d.h"
#include "gassim2d_util.h"

#include <string>
#include <sstream>

#include <iostream>

#include <stdio.h>

#include "simblock.h"
#include "datagen.h"
#include "solution.h"

#include <gtest/gtest.h>
#include "gtest-mpi-listener.hpp"
#include <mpi.h>


#include <string>
#include <sstream>


#include "gas_constants.h"

#include "mpi_gassim_fixture.h"

#define DEFAULT_GENERATOR_SEED 1337

#define TESTING_TOLERANCE_GHOSTS 1E-12
#define TESTING_TOLERANCE 1E-8

bool particle_sort_cmp(phys_particle_t a, phys_particle_t b){
	if(a.p.x < b.p.x) return true;
	if(a.p.x > b.p.x) return false;
	if(a.p.y < b.p.y) return true;
	if(a.p.y > b.p.y) return false;
	if(a.v.x < b.v.x) return true;
	if(a.v.x > b.v.x) return false;
	//if(a.v.y < b.v.y) return true;
	//if(a.v.y > b.v.y) return false;
	return (a.v.y < b.v.y);
}

TEST_F(MpiGasSimBounded, basic_init){

	//fprintf(stderr,"Hello from %i (%i): %i \n", testing_rank, testing_rank, in_testing_group);
	//ASSERT_EQ(true,true);

	if(this->in_testing_group){
		std::shared_ptr<SimulationBlock> my_block = this->create_sim_block();
		my_block->init_communication();
		my_block->finalize_communication();
	}
}

TEST_F(MpiGasSimBounded, part1_GhostsCommunicated_1){
	if(!in_testing_group){
		return;
	}

	const int N_ghost_tests = 8;
	const int expect_ghost_number[][4] = {{0,0,0,0},
									{0,0,0,0},
									{0,1,1,1},
									{1,0,1,1},
									{1,1,0,1},
									{1,1,1,0},
									{0,0,0,0},
									{0,0,0,0}
								};
	//const int ghost_sender[] = {0,1,0,1,2,3,2,3};
	const int ghost_sender[] = {0,1,0,1,2,3,2,3};//-1 will skip that subtest
	const phys_particle_t test_ghosts[] = {
		{{0.1, 0.1}, {0.0, 0.0}},//First ghost, NW, tests both negative wraparounds
		{{19.9, 0.1}, {0.0, 0.0}},//Second ghost, NE, tests positive and negative wraparounds
		{{9.9, 9.9}, {0.0, 0.0}},//Tests SE without wraparound
		{{10.1, 9.9}, {0.0, 0.0}},//Tests SW without wraparound
		{{9.9, 10.1}, {0.0, 0.0}},//Tests NE without wraparound
		{{10.1, 10.1}, {0.0, 0.0}},//Tests NE without wraparound
		{{0.1, 19.9}, {0.0, 0.0}},//Tests SW with positive and negative wraparound
		{{19.9, 19.9}, {0.0, 0.0}},//Tests SE with both positive wraparound
	};

	const phys_vector_t expected_ghost_positions[][4] = {
		{ {0.1, 0.1}, {20.1, 0.1}, {0.1, 20.1}, {20.1, 20.1} },
		{ {-0.1, 0.1}, {19.9, 0.1}, {-0.1, 20.1}, {19.9,20.1} },
		{ {9.9, 9.9}, {9.9, 9.9}, {9.9, 9.9}, {9.9, 9.9} },
		{ {10.1, 9.9}, {10.1, 9.9}, {10.1, 9.9}, {10.1, 9.9} },
		{ {9.9, 10.1}, {9.9, 10.1}, {9.9, 10.1}, {9.9, 10.1} },
		{ {10.1, 10.1}, {10.1, 10.1}, {10.1, 10.1}, {10.1, 10.1} },
		{ {0.1, -0.1}, {20.1, -0.1}, {0.1, 19.9}, {20.1, 19.9} },
		{ {-0.1, -0.1}, {19.9, -0.1}, {-0.1, 19.9}, {19.9, 19.9} },
	};

	std::shared_ptr<SimulationBlock> my_sim_block = this->create_sim_block();
	my_sim_block->init_communication();

	int k;
	for(k=0;k<N_ghost_tests;k++){
		if(-1 == ghost_sender[k]){ continue; }
		MPI_Barrier(testing_comm);
		//clear all particles
		while(my_sim_block->N_particles > 0){ my_sim_block->remove_particle(0); }
		my_sim_block->N_ghosts = 0;
		ASSERT_EQ(0, my_sim_block->N_particles) << "No simblock should have any particles at this point.";
		ASSERT_EQ(0, my_sim_block->N_ghosts) << "No simblock should have any ghosts at this point.";


		if(testing_rank == ghost_sender[k]){
			my_sim_block->add_particle(test_ghosts[k]);
			EXPECT_EQ(1, my_sim_block->N_particles) << "The originating simulation block should have 1 partical";
		}
		MPI_Barrier(testing_comm);
		my_sim_block->communicate_ghosts();
		MPI_Barrier(testing_comm);

		EXPECT_EQ(expect_ghost_number[k][testing_rank],my_sim_block->N_ghosts);

		if(testing_rank != ghost_sender[k] && expect_ghost_number[k][testing_rank] > 0 && my_sim_block->N_ghosts > 0){
			phys_vector_t an_expected_ghost_position = expected_ghost_positions[k][testing_rank];
			phys_particle_t one_ghost = my_sim_block->all_ghosts[0];
			EXPECT_NEAR(an_expected_ghost_position.x, one_ghost.p.x,TESTING_TOLERANCE_GHOSTS);
			EXPECT_NEAR(an_expected_ghost_position.y, one_ghost.p.y,TESTING_TOLERANCE_GHOSTS);
		}
	}



	my_sim_block->finalize_communication();
	my_sim_block.reset();

}

TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_VELOCITY){
	if(!in_testing_group){
		return;
	}

	const phys_particle_t one_test_particle = {.p = {.x = 10.1, .y = 5.0}, .v = {.x = 1234.5, .y = 3495.0}};

	std::shared_ptr<SimulationBlock> my_sim_block = this->create_sim_block();
	my_sim_block->init_communication();

	//clear all particles
	while(my_sim_block->N_particles > 0){ my_sim_block->remove_particle(0); }
	my_sim_block->N_ghosts = 0;
	ASSERT_EQ(0, my_sim_block->N_particles) << "No simblock should have any particles at this point.";
	ASSERT_EQ(0, my_sim_block->N_ghosts) << "No simblock should have any ghosts at this point.";

	if(testing_rank == 0){
		my_sim_block->add_particle(one_test_particle);
	}
	my_sim_block->exchange_particles();

	if(testing_rank == 1){
		EXPECT_EQ(1, my_sim_block->N_particles) << "Did not migrate the correct number of particles";
		if(my_sim_block->N_particles > 0){
			const phys_particle_t received = my_sim_block->all_particles[0];
			EXPECT_EQ(one_test_particle.v.x, received.v.x) << "Did not send correct x velocity when migrate single particle.";
			EXPECT_EQ(one_test_particle.v.y, received.v.y) << "Did not send correct y velocity when migrate single particle.";
		}
	}

	if(testing_rank == 0 ){
		EXPECT_EQ(0, my_sim_block->N_particles) << "Did not remove particles when migrating.";
	}

	if(testing_rank == 2 || testing_rank == 3){
		EXPECT_EQ(0, my_sim_block->N_particles) << "Sent a migrant particle to the wrong destination.";
	}
	EXPECT_EQ(0, my_sim_block->N_ghosts) << "No simblock should have any ghosts at this point.";

	my_sim_block->finalize_communication();

	MPI_Barrier(testing_comm);
	my_sim_block.reset();
	//EXPECT_EQ(true,true);

}

TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_SELF){
	//exchange_tester(4);
	exchange_test_helper(0,0,{5.0,5.0},{5.0,5.0});
}

TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_NW){
	//exchange_tester(12);
	exchange_test_helper(3,0,{9.9,9.9},{9.9,9.9});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_N){
	//exchange_tester(13);
	exchange_test_helper(2,0,{5.0,9.9},{5.0,9.9});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_NE){
	//exchange_tester(14);
	exchange_test_helper(2,1,{10.1,9.9},{10.1,9.9});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_W){
	//exchange_tester(15);
	exchange_test_helper(1,0,{9.9,5.0},{9.9,5.0});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_E){
	//exchange_tester(5);{10.1,5.0}
	exchange_test_helper(0,1,{10.1,5.0},{10.1,5.0});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_SW){
	//exchange_tester(16);
	exchange_test_helper(1,2,{9.9,10.1},{9.9,10.1});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_S){
	//exchange_tester(7);
	exchange_test_helper(0,2,{5.0,10.1},{5.0,10.1});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_SE){
	//exchange_tester(8);
	exchange_test_helper(0,3,{10.1,10.1},{10.1,10.1});
}

TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_NW_boundary){
	//exchange_tester(0);
	boundary_test_helper(0,{-0.1,-0.1});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_N_boundary){
	//exchange_tester(1);
	boundary_test_helper(0,{5.0,-0.1});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_NE_boundary){
	//exchange_tester(2);
	boundary_test_helper(1,{20.1,-0.1});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_W_boundary){
	//exchange_tester(3);
	boundary_test_helper(0,{-0.1,5.0});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_E_boundary){
	//exchange_tester(9);
	boundary_test_helper(1,{20.1,5.0});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_SW_boundary){
	//exchange_tester(6);
	boundary_test_helper(2,{-0.1,20.1});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_S_boundary){
	//exchange_tester(10);
	boundary_test_helper(2,{5.0,20.1});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_SE_boundary){
	//exchange_tester(11);
	boundary_test_helper(3,{20.1,20.1});
}

TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_NW_slide_W){
	//exchange_tester(0);
	exchange_test_helper(1,0,{9.9,-0.1},{9.9,-0.1});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_NW_slide_N){
	//exchange_tester(0);
	exchange_test_helper(2,0,{-0.1,9.9},{-0.1,9.9});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_NE_slide_E){
	//exchange_tester(2);
	exchange_test_helper(0,1,{10.1,-0.1},{10.1,-0.1});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_NE_slide_N){
	//exchange_tester(2);
	exchange_test_helper(3,1,{20.1,9.9},{20.1,9.9});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_SW_slide_W){
	//exchange_tester(6);
	exchange_test_helper(3,2,{9.9,20.1},{9.9,20.1});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_SW_slide_S){
	//exchange_tester(6);
	exchange_test_helper(0,2,{-0.1,10.1},{-0.1,10.1});
}

TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_SE_slide_E){
	//exchange_tester(6);
	exchange_test_helper(2,3,{10.1,20.1},{10.1,20.1});
}
TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_SE_slide_S){
	//exchange_tester(6);
	exchange_test_helper(1,3,{20.1,10.1},{20.1,10.1});
}


TEST_F(MpiGasSimBounded, part1_ParticlesExchanged_MULTIPLE){
	if(!in_testing_group){
		return;
	}


	const int N_exchange_multiple = 3;
	const phys_particle_t test_particles[] = {
		{{10.1, 5.0}, {12.30, -345.6}},
		{{10.11, 5.0}, {-777.0, 292.0}},
		{{10.2, 5.0}, {88.0, -90.1}},
	};



	std::shared_ptr<SimulationBlock> my_sim_block = this->create_sim_block();
	my_sim_block->init_communication();

	//clear all particles
	while(my_sim_block->N_particles > 0){ my_sim_block->remove_particle(0); }
	my_sim_block->N_ghosts = 0;
	ASSERT_EQ(0, my_sim_block->N_particles) << "No simblock should have any particles at this point.";
	ASSERT_EQ(0, my_sim_block->N_ghosts) << "No simblock should have any ghosts at this point.";

	if(testing_rank == 0){
		for(int j = 0;j<N_exchange_multiple;++j){
			my_sim_block->add_particle(test_particles[j]);
		}
	}
	my_sim_block->exchange_particles();

	if(testing_rank == 1){
		EXPECT_EQ(N_exchange_multiple, my_sim_block->N_particles) << "Did not communicate the correct number of ghosts";
		//Student may not have communicated them in the order I had them in.
		//That is perfectly acceptable. But we need to confirm they are all there.
		std::sort(my_sim_block->all_particles, my_sim_block->all_particles + my_sim_block->N_particles, particle_sort_cmp);

		for(int j = 0;j<N_exchange_multiple;++j){
			phys_particle_t expected_particle = test_particles[j];
			phys_particle_t migrant = my_sim_block->all_particles[j];

			EXPECT_NEAR(migrant.p.x, expected_particle.p.x, TESTING_TOLERANCE_GHOSTS) << "Migrant .p.x wrong.";
			EXPECT_NEAR(migrant.p.y, expected_particle.p.y, TESTING_TOLERANCE_GHOSTS) << "Migrant .p.y wrong.";
			EXPECT_NEAR(migrant.v.x, expected_particle.v.x, TESTING_TOLERANCE_GHOSTS) << "Migrant .v.x wrong.";
			EXPECT_NEAR(migrant.v.y, expected_particle.v.y, TESTING_TOLERANCE_GHOSTS) << "Migrant .v.y wrong.";
		}
	}

	if(testing_rank == 0 || testing_rank == 2 || testing_rank == 3){
		EXPECT_EQ(0, my_sim_block->N_particles) << "Did not communicate the correct number of particles";
	}

	my_sim_block->finalize_communication();

	MPI_Barrier(testing_comm);
	my_sim_block.reset();
	//EXPECT_EQ(true,true);

}

/*
TEST_F(MpiGasSimBounded, part1_ParticlesMotion){
	if(!in_testing_group){
		return;
	}

	std::shared_ptr<SimulationBlock> my_sim_block = this->create_sim_block();

	my_sim_block->tick_delta_t_ns = 0.01;//No collisions or interactions so we can have bigger ticks.
	int num_steps = (int)( 4.0/my_sim_block->tick_delta_t_ns );
	my_sim_block->exchange_frequency = (int)(num_steps/8);
	if(my_sim_block->exchange_frequency < 1) my_sim_block->exchange_frequency = 1;

	if(testing_rank == 0){
		phys_particle_t par;
		par.p.x = 0.0625;
		par.p.y = 0.0625;
		par.v.x = 20.0;
		par.v.y = 10.0;

		my_sim_block->add_particle(par);
	}

	my_sim_block->init_communication();

	int i;
	for(i=0;i<num_steps;i++){
		my_sim_block->simulate_tick();
	}
	my_sim_block->exchange_particles();//Ensure that a final particle echange happens.

	if(testing_rank == 0){
		EXPECT_EQ(1,my_sim_block->N_particles);

		EXPECT_NEAR(0.0625,my_sim_block->all_particles[0].p.x,TESTING_TOLERANCE);
		EXPECT_NEAR(0.0625,my_sim_block->all_particles[0].p.y,TESTING_TOLERANCE);

	}else{
		EXPECT_EQ(0,my_sim_block->N_particles);
	}

	//
	//Tests migration in other directions
	//

	if(testing_rank == 0){
		my_sim_block->remove_particle(0);

		phys_particle_t par;
		par.p.x = 0.0625;
		par.p.y = 0.0625;
		par.v.x = -20.0;
		par.v.y = -10.0;

		my_sim_block->add_particle(par);
	}

	//my_sim_block->init_communication();

	for(i=0;i<num_steps;i++){
		my_sim_block->simulate_tick();
	}
	my_sim_block->exchange_particles();//Ensure that a final particle echange happens.

	if(testing_rank == 0){
		EXPECT_EQ(1,my_sim_block->N_particles);

		EXPECT_NEAR(0.0625,my_sim_block->all_particles[0].p.x,TESTING_TOLERANCE);
		EXPECT_NEAR(0.0625,my_sim_block->all_particles[0].p.y,TESTING_TOLERANCE);

	}else{
		EXPECT_EQ(0,my_sim_block->N_particles);
	}

	//
 	//DONE
	//

	my_sim_block->finalize_communication();
	my_sim_block.reset();

}
*/

int main(int argc, char** argv) {
  // Filter out Google Test arguments
  ::testing::InitGoogleTest(&argc, argv);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int world_rank;
  int world_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  if(4 > world_size){
	  std::cerr << " Can't run these tests on fewer than 4 ranks " << std::endl;
	  exit(1);
  }

  // Add object that will finalize MPI on exit; Google Test owns this pointer
  ::testing::AddGlobalTestEnvironment(new MPIEnvironment);

  // Get the event listener list.
  ::testing::TestEventListeners& listeners =
      ::testing::UnitTest::GetInstance()->listeners();

  //only get json/xml output from rank 0
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(0 != rank){
	  // Remove default listener
	  delete listeners.Release(listeners.default_result_printer());
	  // Remove default file output
	  delete listeners.Release(listeners.default_xml_generator());
  }
  listeners.Append(new MPICollectTestResults);

  // Run tests, then clean up and exit
  int foo = RUN_ALL_TESTS();

  return foo;
}
