#include <stdint.h>
#include <string.h>
#include <stdbool.h>

#include "afl-fuzz.h"
#include "types.h"

#ifndef DARWIN_H
#define DARWIN_H

typedef unsigned int uint;
#define REAL_VALUED
#ifdef REAL_VALUED
	#define VALUE_TYPE double
#else
	#define VALUE_TYPE bool
#endif

/*typedef struct {
	unsigned int *_nextToEvaluate;  // next child to evaluate
	unsigned int *_currentParent;   // current parent being mutated
	unsigned int *_bestSoFar;       // best so far parent
	unsigned int *_bestSoFarChild;  // best so far child (for a single parent)
	int **_parentsFitness;
	int **_childrenFitness;

	VALUE_TYPE ***_parents;
	VALUE_TYPE ***_children;
	VALUE_TYPE **_currentRelativeProb;
}cluster_data_structures_for_AMBES;*/

#endif


#ifndef DARWIN_F
#define DARWIN_F
void mesfuzz_DAR_init(afl_state_t *afl, u32 nr_seeds, unsigned nr_mutations, u32 group_id,cluster_data_structures_for_AMBES cluster_evo_data_x[]);

int mesfuzz_DAR_SelectOperator(u32 seed, u32 group_id,cluster_data_structures_for_AMBES cluster_evo_data_x[]);

void mesfuzz_DAR_NotifyFeedback(afl_state_t *afl, u32 seed, unsigned numPaths, u32 group_id,cluster_data_structures_for_AMBES cluster_evo_data_x[]);

u32 mesfuzz_DAR_get_parent_repr(u32 seed, u32 group_id, cluster_data_structures_for_AMBES cluster_evo_data_x[]);



#endif
