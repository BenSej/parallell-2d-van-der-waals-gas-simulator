#include "solution.h"

#include <cassert>
#include <stdexcept>
#include <cstdlib>
#include <iostream>

#include <mpi.h>
#include <list>
#include <vector>




//Do not alter the initializer line, nor the signature of the constructor.
MPISimulationBlock::MPISimulationBlock(gas_simulation _in_sim, SimulationGrid _grid, MPI_Comm _comm, int _rank) : SimulationBlock(_in_sim, _grid, _rank, 10000), my_comm(_comm), my_rank(_rank){
	/* You may add code at the bottom of this constructor */
}

MPISimulationBlock::~MPISimulationBlock(){
}

/* Students overwrite and implement these. */
int MPISimulationBlock::init_communication(){
	return 0;
}
int MPISimulationBlock::finalize_communication(){
	return 0;
}
int MPISimulationBlock::exchange_particles(){

	//NOTE: because "remove_particle" creates a hole and
	//fills it with the last particle in the array, you MUST iterate backwards.
	//(In fact, you will find in general in programming that deletion should
	// usually be done with backwards iterators.
	// I have found bugs in widely-distributed, production NASA software
	// that arise from deleting on a forward iterator.)
	//     Loop signature should be similar to:
	//
	//    for(int i=N_particles-1;i>=0;--i){
	//		FIGURE OUT PARTICLE MIGRATION DIRECTION
	//	}
	//
	//	MIGRATE ALL PARTICLES

	SimulationGrid::grid_position my_pos = the_grid.get_grid_position(block_rank);

	std::vector<phys_particle_t> particles[8];
	MPI_Request recv_requests_num_particles[8];
	MPI_Request send_requests_num_particles[8];
	MPI_Request recv_requests_particles[8];
	MPI_Request send_requests_particles[8];

	for (int i = N_particles - 1; i >= 0; i--) {
	    phys_particle_t p = all_particles[i];
	    sim_direction_t direction = check_migrant_direction(p);
		if (direction == DIR_SELF) continue;

		switch(direction) {
			case DIR_N:
				particles[0].push_back(p);
				break;
			case DIR_S:
				particles[1].push_back(p);
				break;
			case DIR_E:
				particles[2].push_back(p);
				break;
			case DIR_W:
				particles[3].push_back(p);
				break;
			case DIR_NE:
				particles[4].push_back(p);
				break;
			case DIR_NW:
				particles[5].push_back(p);
				break;
			case DIR_SE:
				particles[6].push_back(p);
				break;
			case DIR_SW:
				particles[7].push_back(p);
				break;
		}
		remove_particle(i);
	}


	for (int i = 0; i < 8; i++) {
		SimulationGrid::grid_position neighbor_position = my_pos;
		switch(i) {
			case 0:
				neighbor_position.i--;
				break;
			case 1:
				neighbor_position.i++;
				break;
			case 2:
				neighbor_position.j++;
				break;
			case 3:
				neighbor_position.j--;
				break;
			case 4:
				neighbor_position.i--;
				neighbor_position.j++;
				break;
			case 5:
				neighbor_position.i--;
				neighbor_position.j--;
				break;
			case 6:
				neighbor_position.i++;
				neighbor_position.j++;
				break;
			case 7:
				neighbor_position.i++;
				neighbor_position.j--;
				break;
		}
		int rank = the_grid.cpu_for_position(neighbor_position.i, neighbor_position.j);

		if (rank == my_rank) continue;

		int num_particles_send = particles[i].size();

		// tell this direction how many particles we're sending over
		MPI_Isend(&num_particles_send, 1, MPI_INT, rank, 0, my_comm, &send_requests_num_particles[i]);

		// send the particles to this direction
		if (num_particles_send > 0) {
			MPI_Isend(particles[i].data(), num_particles_send * sizeof(phys_particle_t), MPI_BYTE, rank, 0, my_comm, &send_requests_particles[i]);
		}

	}

	for (int i = 0; i < 8; i++) {
		SimulationGrid::grid_position neighbor_position = my_pos;
		switch(i) {
			case 0:
				neighbor_position.i--;
				break;
			case 1:
				neighbor_position.i++;
				break;
			case 2:
				neighbor_position.j++;
				break;
			case 3:
				neighbor_position.j--;
				break;
			case 4:
				neighbor_position.i--;
				neighbor_position.j++;
				break;
			case 5:
				neighbor_position.i--;
				neighbor_position.j--;
				break;
			case 6:
				neighbor_position.i++;
				neighbor_position.j++;
				break;
			case 7:
				neighbor_position.i++;
				neighbor_position.j--;
				break;
		}
		int rank = the_grid.cpu_for_position(neighbor_position.i, neighbor_position.j);

		if (rank == my_rank) continue;

		int num_particles_recv = 0;
		// receive the message from this direction telling us how many particles to receive
		MPI_Irecv(&num_particles_recv, 1, MPI_INT, rank, 0, my_comm, &recv_requests_num_particles[i]);
		MPI_Wait(&recv_requests_num_particles[i], MPI_STATUS_IGNORE);

		// receive the particles from this direction'
		std::vector<phys_particle_t> particles_recv(num_particles_recv);
		if (num_particles_recv > 0) {
			MPI_Irecv(particles_recv.data(), num_particles_recv * sizeof(phys_particle_t), MPI_BYTE, rank, 0, my_comm, &recv_requests_particles[i]);
			MPI_Wait(&recv_requests_particles[i], MPI_STATUS_IGNORE);
		}

		// add particles we received
		for (phys_particle_t p : particles_recv) {
			add_particle(p);
		}
	}

	return 0;
}

int MPISimulationBlock::communicate_ghosts(){
	/*
	* Our superclass (defined in "simblock.h")
	* has (and therefore we have) the member variables:
	*
	* `unsigned int N_ghosts`
	* and
	* `phys_particle_t *all_ghosts`
	*
	* In this function, you should find all the particles that need to be
	* communicated as ghosts, and communicate them (or at least their position)
	* to adjacent blocks.
	*
	* The particles that you recieve from adjacent blocks need to be added
	* to `this->all_ghosts` and `this->N_ghosts` needs to be changed to equal
	* the total number of ghosts received in this call.
	*
	* Ghosts are reset every time this is called, so you should start populating
	* `this->all_ghosts` from the 0th position.
	*
	* If this block is on a Northern, Southern, Eastern, or Western edge,
	* make sure to send to the appropriate wrapped block. In this case,
	* You will also need to adjust the X or Y coordinate
	* (`p.x`, or `p.y` members of the `phys_particle_t` type.)
	* accordingly for wraparound.
	*
	*/

	this->N_ghosts = 0;//Clear all the ghosts we have, we get entirely new ghosts.
	//TODO: Add your solution after this line.

	SimulationGrid::grid_position my_pos = the_grid.get_grid_position(block_rank);

	std::vector<phys_particle_t> particles[8];
	MPI_Request recv_requests_num_particles[8];
	MPI_Request send_requests_num_particles[8];
	MPI_Request recv_requests_particles[8];
	MPI_Request send_requests_particles[8];

	for (int i = N_particles - 1; i >= 0; i--) {
	    phys_particle_t p = all_particles[i];
	    sim_direction_t direction = check_ghost_direction(p);
		if (direction == DIR_SELF) continue;

		if (DIR_HAS(DIR_N, direction)) {
			particles[0].push_back(p);
		}
		if (DIR_HAS(DIR_S, direction)) {
			particles[1].push_back(p);
		}
		if (DIR_HAS(DIR_E, direction)) {
			particles[2].push_back(p);
		}
		if (DIR_HAS(DIR_W, direction)) {
			particles[3].push_back(p);
		}
		if (DIR_HAS(DIR_NE, direction)) {
			particles[4].push_back(p);
		}
		if (DIR_HAS(DIR_NW, direction)) {
			particles[5].push_back(p);
		}
		if (DIR_HAS(DIR_SE, direction)) {
			particles[6].push_back(p);
		}
		if (DIR_HAS(DIR_SW, direction)) {
			particles[7].push_back(p);
		}
	}

	for (int i = 0; i < 8; i++) {
		SimulationGrid::grid_position neighbor_position = my_pos;
		switch(i) {
			case 0:
				neighbor_position.i++;
				break;
			case 1:
				neighbor_position.i--;
				break;
			case 2:
				neighbor_position.j++;
				break;
			case 3:
				neighbor_position.j--;
				break;
			case 4:
				neighbor_position.i++;
				neighbor_position.j++;
				break;
			case 5:
				neighbor_position.i++;
				neighbor_position.j--;
				break;
			case 6:
				neighbor_position.i--;
				neighbor_position.j++;
				break;
			case 7:
				neighbor_position.i--;
				neighbor_position.j--;
				break;
		}

		int rank = the_grid.cpu_for_position(neighbor_position.i, neighbor_position.j);
		if (rank == my_rank) continue;

		int num_particles_send = particles[i].size();

		// tell this direction how many particles we're sending over
		MPI_Isend(&num_particles_send, 1, MPI_INT, rank, 0, my_comm, &send_requests_num_particles[i]);

		// send the particles to this direction
		if (num_particles_send > 0) {
			MPI_Isend(particles[i].data(), num_particles_send * sizeof(phys_particle_t), MPI_BYTE, rank, 0, my_comm, &send_requests_particles[i]);
		}

	}
	


	for (int i = 0; i < 8; i++) {
		SimulationGrid::grid_position neighbor_position = my_pos;
		switch(i) {
			case 0:
				neighbor_position.i++;
				break;
			case 1:
				neighbor_position.i--;
				break;
			case 2:
				neighbor_position.j++;
				break;
			case 3:
				neighbor_position.j--;
				break;
			case 4:
				neighbor_position.i++;
				neighbor_position.j++;
				break;
			case 5:
				neighbor_position.i++;
				neighbor_position.j--;
				break;
			case 6:
				neighbor_position.i--;
				neighbor_position.j++;
				break;
			case 7:
				neighbor_position.i--;
				neighbor_position.j--;
				break;
		}

		int rank = the_grid.cpu_for_position(neighbor_position.i, neighbor_position.j);
		if (rank == my_rank) continue;

		int num_particles_recv = 0;
		// receive the message from this direction telling us how many particles to receive
		MPI_Irecv(&num_particles_recv, 1, MPI_INT, rank, 0, my_comm, &recv_requests_num_particles[i]);
		MPI_Wait(&recv_requests_num_particles[i], MPI_STATUS_IGNORE);

		// receive the particles from this direction
		std::vector<phys_particle_t> particles_recv(num_particles_recv);
		if (num_particles_recv > 0) {
			MPI_Irecv(particles_recv.data(), num_particles_recv * sizeof(phys_particle_t), MPI_BYTE, rank, 0, my_comm, &recv_requests_particles[i]);
			MPI_Wait(&recv_requests_particles[i], MPI_STATUS_IGNORE);
		}

		// add particles we received
		for (int j = 0; j < particles_recv.size(); j++) {
			this->all_ghosts[j + this->N_ghosts] = particles_recv[j];
		}
		this->N_ghosts += num_particles_recv;
	}


	return 0;
}
