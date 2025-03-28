#ifndef P1_SOLLUTION_H
#define P1_SOLLUTION_H

#include "gassim2d.h"
#include "simblock.h"

#include <mpi.h>
//Include whatever you like.

//You may define your own helper classes if you wish.

class MPISimulationBlock : public SimulationBlock {
public:
	//Don't alter the signature of the constructor.
	MPISimulationBlock(gas_simulation _in_sim, SimulationGrid _grid, MPI_Comm _comm, int _rank);
	~MPISimulationBlock();//Make sure to clean up after yourself.

	/* Students overwrite and implement these. */
	virtual int init_communication();
	virtual int finalize_communication();
	virtual int exchange_particles();
	virtual int communicate_ghosts();
	
private://but students can access
	//SimulationGrid the_grid;
	MPI_Comm my_comm;
	int my_rank;

	//You may need to add your own variables.
};





#endif /* P1_SOLLUTION_H */
