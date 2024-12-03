
#include "afl-darwin.h"
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "afl-fuzz.h"


unsigned kAflNumMutationOperators = 22;

#define REAL_VALUED
#ifdef REAL_VALUED
	#define VALUE_TYPE double
#else
	#define VALUE_TYPE bool
#endif

#ifndef _nParents
#define _nParents 5u
#endif

#ifndef _lambda
#define _lambda 4u
#endif


/* The largest number rand will return (same as INT_MAX).  */
#define	RAND_MAX	2147483647

void rand_init();
uint64_t romuDuoJr_random ();
double rand_32_double();
unsigned int rand_32_int (unsigned int max);
double rand_32_double_gauss();
double generate_uniform_random();

/*
 * Initialize data structures for 
 * @nr_seeds: number of different initial seeds
 * @nr_mutations: number of mutation operators
*/

void mesfuzz_DAR_init(afl_state_t *afl, u32 nr_seeds, unsigned nr_mutations, u32 group_id,cluster_data_structures_for_AMBES cluster_evo_data_x[]) {
	// init RNG
	rand_init();

	nr_seeds = 1;
	group_id -= 1;

	kAflNumMutationOperators = nr_mutations;

	// initialize opt alg data structures
	
	cluster_evo_data_x[group_id]._nextToEvaluate = (unsigned int *) malloc(nr_seeds * sizeof(unsigned int));
	cluster_evo_data_x[group_id]._currentParent = (unsigned int *) malloc(nr_seeds * sizeof(unsigned int));
	cluster_evo_data_x[group_id]._bestSoFar = (unsigned int *) malloc(nr_seeds * sizeof(unsigned int));         // best so far parent
	cluster_evo_data_x[group_id]._bestSoFarChild = (unsigned int *) malloc(nr_seeds * sizeof(unsigned int));    // best so far child (for a single parent)
	cluster_evo_data_x[group_id]._parentsFitness = (int **) malloc(nr_seeds * sizeof(int *));
	cluster_evo_data_x[group_id]._childrenFitness = (int **) malloc(nr_seeds * sizeof(int *));

	cluster_evo_data_x[group_id]._parents = (VALUE_TYPE ***) malloc(nr_seeds * sizeof(VALUE_TYPE **));
	cluster_evo_data_x[group_id]._children = (VALUE_TYPE ***) malloc(nr_seeds * sizeof(VALUE_TYPE **));
	
	cluster_evo_data_x[group_id]._currentRelativeProb = (VALUE_TYPE **) malloc(nr_seeds * sizeof(VALUE_TYPE *));
	

	//printf("nParents: %u, lambda %u\n", _nParents, _lambda);

	for (int seed = 0; seed < nr_seeds; seed++) {
		cluster_evo_data_x[group_id]._nextToEvaluate[seed] = 0;
		cluster_evo_data_x[group_id]._bestSoFar[seed] = 0;
		cluster_evo_data_x[group_id]._bestSoFarChild[seed] = 0;
		cluster_evo_data_x[group_id]._currentParent[seed] = 0;

		cluster_evo_data_x[group_id]._parentsFitness[seed] = malloc(_nParents * sizeof(int));
		cluster_evo_data_x[group_id]._childrenFitness[seed] = malloc(_lambda * sizeof(int));

		cluster_evo_data_x[group_id]._children[seed] = (VALUE_TYPE **) malloc(_lambda * sizeof(VALUE_TYPE *));
		for (unsigned int i = 0; i < _lambda; i++)
			cluster_evo_data_x[group_id]._children[seed][i] = malloc(kAflNumMutationOperators * sizeof(VALUE_TYPE));

		cluster_evo_data_x[group_id]._parents[seed] = (VALUE_TYPE **) malloc(_nParents * sizeof(VALUE_TYPE *));
		for (unsigned int i = 0; i < _nParents; i++)
			cluster_evo_data_x[group_id]._parents[seed][i] = malloc(kAflNumMutationOperators * sizeof(VALUE_TYPE));

		memset(cluster_evo_data_x[group_id]._parentsFitness[seed], 0, _nParents);
		memset(cluster_evo_data_x[group_id]._childrenFitness[seed], 0, _lambda);

		// initial random values for the parents and first child individual
		for(unsigned int i = 0; i < _nParents; i++) {
			for(unsigned int j = 0; j < kAflNumMutationOperators; j++) {
				// TODO: optionally replace with an alternate randomizer
				//cluster_evo_data_x[group_id]._parents[seed][i][j] = rand() > (RAND_MAX / 2);
				cluster_evo_data_x[group_id]._parents[seed][i][j] = 0.5;
				unsigned int mm = cluster_evo_data_x[group_id]._parents[seed][i][j];
			}
		}
		double *v = cluster_evo_data_x[group_id]._parents[seed][0];
		creat_cluster_prob_file(afl, 0, group_id,cluster_evo_data_x,22);
	
		for(unsigned int i = 0; i < kAflNumMutationOperators; i++) {
			cluster_evo_data_x[group_id]._children[seed][cluster_evo_data_x[group_id]._nextToEvaluate[seed]][i] = cluster_evo_data_x[group_id]._parents[seed][0][i] + (rand_32_double_gauss() - 0.5) / 4 *10;
		}
		int randomGene = rand_32_int(kAflNumMutationOperators);
		double mutation = (rand_32_double_gauss() - 0.5) / 4 * 10;
		cluster_evo_data_x[group_id]._children[seed][cluster_evo_data_x[group_id]._nextToEvaluate[seed]][randomGene] += mutation;
		if(cluster_evo_data_x[group_id]._children[seed][cluster_evo_data_x[group_id]._nextToEvaluate[seed]][randomGene] < 0)
			cluster_evo_data_x[group_id]._children[seed][cluster_evo_data_x[group_id]._nextToEvaluate[seed]][randomGene] = 0;

		double *first_child = cluster_evo_data_x[group_id]._children[seed][cluster_evo_data_x[group_id]._nextToEvaluate[seed]];

		cluster_evo_data_x[group_id]._currentRelativeProb[seed] = cluster_evo_data_x[group_id]._children[seed][cluster_evo_data_x[group_id]._nextToEvaluate[seed]];
	}
}

/*
 * Choose an AFL mutation operator
 * @seed: seed to select per-seed vector
*/
int MesFuzz_SelectOperator(u32 seed, u32 group_id,cluster_data_structures_for_AMBES cluster_evo_data_x[])
{
	group_id -= 1;
#ifdef REAL_VALUED
	cluster_data_structures_for_AMBES cluster_date = cluster_evo_data_x[group_id];
	double *v = cluster_evo_data_x[group_id]._currentRelativeProb[seed];
	uint v_size = kAflNumMutationOperators;
	double cumulative[kAflNumMutationOperators];
	for (uint i = 0; i < kAflNumMutationOperators; i++)	
		cumulative[i] = 0.0;

	double minVal = 0.0;
	double maxVal = 1.0;
	for (uint i = 0; i < kAflNumMutationOperators; i++) {
		if (v[kAflNumMutationOperators] < minVal) {
			minVal = v[kAflNumMutationOperators];
		} else {
			if (v[kAflNumMutationOperators] > maxVal) {
				maxVal = v[kAflNumMutationOperators];
			}
		}
	}

	cumulative[0] = 1 + (kAflNumMutationOperators - 1) * (v[0] - minVal) / (maxVal - minVal); // selection pressure is kAflNumMutationOperators
	for(uint i = 1; i < kAflNumMutationOperators; i++) {
		cumulative[i] = 1 + (kAflNumMutationOperators - 1) * (v[i] - minVal)/(maxVal - minVal);
		cumulative[i] += cumulative[i - 1];
	}

	double rngNormVal = rand_32_double();

	double randVal = rngNormVal * cumulative[kAflNumMutationOperators - 1];

	int chosen = 0;
	// didn't bother with binary search since negligible compared to other modules
	while(cumulative[chosen] < randVal)
		chosen++;
	return chosen;
#else
	// baseline:
	bool *v = (bool *) cluster_evo_data_x[group_id]._currentRelativeProb[seed];

	// select a random mutation operator with flag == true
	// revert to random if all false
	unsigned int operatorId = rand_32_int(kAflNumMutationOperators);
	unsigned int nTries = 0;
	while (v[operatorId] == false && nTries < kAflNumMutationOperators) {
		nTries++;
		operatorId = (operatorId + 1) % kAflNumMutationOperators;
	}
    	return operatorId;
#endif
}


/*
 * Report feedback to
 * @seed: seed to attribute this to
 * @numPaths: number of new paths
*/
void MesFuzz_NotifyFeedback(afl_state_t *afl, u32 seed, unsigned numPaths, u32 group_id,cluster_data_structures_for_AMBES cluster_evo_data_x[])
{
	group_id -= 1;
	// update this candidate solution fitness
	cluster_evo_data_x[group_id]._childrenFitness[seed][cluster_evo_data_x[group_id]._nextToEvaluate[seed]] = numPaths;

	// update if new best child found
	if(cluster_evo_data_x[group_id]._childrenFitness[seed][cluster_evo_data_x[group_id]._nextToEvaluate[seed]] > cluster_evo_data_x[group_id]._childrenFitness[seed][cluster_evo_data_x[group_id]._bestSoFarChild[seed]]) {
		cluster_evo_data_x[group_id]._bestSoFarChild[seed] = cluster_evo_data_x[group_id]._nextToEvaluate[seed];
	}
	// move to next child candidate
	cluster_evo_data_x[group_id]._nextToEvaluate[seed]++;
	
	// if all children evaluated
	if(cluster_evo_data_x[group_id]._nextToEvaluate[seed] == _lambda) {
		// set best child as future parent (only if better than the parent)
		if(cluster_evo_data_x[group_id]._childrenFitness[seed][cluster_evo_data_x[group_id]._bestSoFarChild[seed]] >cluster_evo_data_x[group_id]._parentsFitness[seed][cluster_evo_data_x[group_id]._currentParent[seed]]) {
			for(unsigned int i = 0; i < kAflNumMutationOperators; i++) {
				cluster_evo_data_x[group_id]._parents[seed][cluster_evo_data_x[group_id]._currentParent[seed]][i] = cluster_evo_data_x[group_id]._children[seed][cluster_evo_data_x[group_id]._bestSoFarChild[seed]][i];
			}
			
			double *v = cluster_evo_data_x[group_id]._parents[seed][cluster_evo_data_x[group_id]._currentParent[seed]];
		  	print_cluster_prob(afl,  0, group_id,cluster_evo_data_x,22);

			cluster_evo_data_x[group_id]._parentsFitness[seed][cluster_evo_data_x[group_id]._currentParent[seed]] = cluster_evo_data_x[group_id]._childrenFitness[seed][cluster_evo_data_x[group_id]._bestSoFarChild[seed]];
		}
		// update best parent solution if needed
		if(cluster_evo_data_x[group_id]._parentsFitness[seed][cluster_evo_data_x[group_id]._currentParent[seed]] > cluster_evo_data_x[group_id]._parentsFitness[seed][cluster_evo_data_x[group_id]._bestSoFar[seed]]) {
			cluster_evo_data_x[group_id]._bestSoFar[seed] = cluster_evo_data_x[group_id]._currentParent[seed];
		}

		// move to next parent (or return to first)
		cluster_evo_data_x[group_id]._currentParent[seed] = (cluster_evo_data_x[group_id]._currentParent[seed] + 1) % _nParents;

		// reset indices
		cluster_evo_data_x[group_id]._bestSoFarChild[seed] = 0;
		cluster_evo_data_x[group_id]._nextToEvaluate[seed] = 0;

		// reset children scores
		memset(cluster_evo_data_x[group_id]._childrenFitness[seed], 0, _lambda);
		
	}
	
	if(cluster_evo_data_x[group_id]._nextToEvaluate[seed] < _lambda) {
		cluster_evo_data_x[group_id]._currentRelativeProb[seed] = &(cluster_evo_data_x[group_id]._children[seed][cluster_evo_data_x[group_id]._nextToEvaluate[seed]][0]);

		kAflNumMutationOperators = 22;
		for(unsigned int i = 1; i < kAflNumMutationOperators; i++) {
			cluster_evo_data_x[group_id]._currentRelativeProb[seed][i] = cluster_evo_data_x[group_id]._parents[seed][cluster_evo_data_x[group_id]._currentParent[seed]][i];
		}
		 double * curr_parent = cluster_evo_data_x[group_id]._parents[seed][cluster_evo_data_x[group_id]._currentParent[seed]];
		
		for(int i =1; i <= 6; i++)
		{
			int randomGene = rand_32_int(kAflNumMutationOperators);
#ifdef REAL_VALUED
			//double mutation = (rand_32_double_gauss() - 0.5) / 4 *10;
			double mutation = (generate_uniform_random() - 0.5) / 4 *10;
			cluster_evo_data_x[group_id]._currentRelativeProb[seed][randomGene] += mutation;
			if(cluster_evo_data_x[group_id]._currentRelativeProb[seed][randomGene] < 0)
				cluster_evo_data_x[group_id]._currentRelativeProb[seed][randomGene] = 0;
		}
		double *New_child = cluster_evo_data_x[group_id]._currentRelativeProb[seed];
		int huahua;
		/*int randomGene = rand_32_int(kAflNumMutationOperators);
#ifdef REAL_VALUED
		double mutation = (rand_32_double_gauss() - 0.5) / 4;
		cluster_evo_data_x[group_id]._currentRelativeProb[seed][randomGene] += mutation;
		if(cluster_evo_data_x[group_id]._currentRelativeProb[seed][randomGene] < 0)
			cluster_evo_data_x[group_id]._currentRelativeProb[seed][randomGene] = 0;*/
#else		
		cluster_evo_data_x[group_id]._currentRelativeProb[seed][randomGene] = !cluster_evo_data_x[group_id]._currentRelativeProb[seed][randomGene];
#endif
	}
}

void MesFuzz_NotifyFeedback_2(afl_state_t *afl, u32 seed, unsigned numPaths, u32 group_id,cluster_data_structures_for_AMBES cluster_evo_data_x[])
{
	group_id -= 1;
	// update this candidate solution fitness
	cluster_evo_data_x[group_id]._childrenFitness[seed][cluster_evo_data_x[group_id]._nextToEvaluate[seed]] = numPaths;

	// update if new best child found
	if(cluster_evo_data_x[group_id]._childrenFitness[seed][cluster_evo_data_x[group_id]._nextToEvaluate[seed]] > cluster_evo_data_x[group_id]._childrenFitness[seed][cluster_evo_data_x[group_id]._bestSoFarChild[seed]]) {
		cluster_evo_data_x[group_id]._bestSoFarChild[seed] = cluster_evo_data_x[group_id]._nextToEvaluate[seed];
	}
	// move to next child candidate
	cluster_evo_data_x[group_id]._nextToEvaluate[seed]++;
	
	// if all children evaluated
	if(cluster_evo_data_x[group_id]._nextToEvaluate[seed] == _lambda) {
		// set best child as future parent (only if better than the parent)
		if(cluster_evo_data_x[group_id]._childrenFitness[seed][cluster_evo_data_x[group_id]._bestSoFarChild[seed]] >= 1) {
			for(unsigned int i = 0; i < kAflNumMutationOperators; i++) {
				cluster_evo_data_x[group_id]._parents[seed][cluster_evo_data_x[group_id]._currentParent[seed]][i] = cluster_evo_data_x[group_id]._children[seed][cluster_evo_data_x[group_id]._bestSoFarChild[seed]][i];
			}
			
			double *v = cluster_evo_data_x[group_id]._parents[seed][cluster_evo_data_x[group_id]._currentParent[seed]];
		  	print_cluster_prob(afl,  0, group_id,cluster_evo_data_x,22);

			cluster_evo_data_x[group_id]._parentsFitness[seed][cluster_evo_data_x[group_id]._currentParent[seed]] = cluster_evo_data_x[group_id]._childrenFitness[seed][cluster_evo_data_x[group_id]._bestSoFarChild[seed]];
		}
		// update best parent solution if needed
		if(cluster_evo_data_x[group_id]._parentsFitness[seed][cluster_evo_data_x[group_id]._currentParent[seed]] > cluster_evo_data_x[group_id]._parentsFitness[seed][cluster_evo_data_x[group_id]._bestSoFar[seed]]) {
			cluster_evo_data_x[group_id]._bestSoFar[seed] = cluster_evo_data_x[group_id]._currentParent[seed];
		}

		// move to next parent (or return to first)
		cluster_evo_data_x[group_id]._currentParent[seed] = (cluster_evo_data_x[group_id]._currentParent[seed] + 1) % _nParents;

		// reset indices
		cluster_evo_data_x[group_id]._bestSoFarChild[seed] = 0;
		cluster_evo_data_x[group_id]._nextToEvaluate[seed] = 0;

		// reset children scores
		memset(cluster_evo_data_x[group_id]._childrenFitness[seed], 0, _lambda);
		
	}
	
	// if there are children to evaluate, generate new candidate and return
	if(cluster_evo_data_x[group_id]._nextToEvaluate[seed] < _lambda) {
		cluster_evo_data_x[group_id]._currentRelativeProb[seed] = &(cluster_evo_data_x[group_id]._children[seed][cluster_evo_data_x[group_id]._nextToEvaluate[seed]][0]);

		kAflNumMutationOperators = 22;
		for(unsigned int i = 0; i < kAflNumMutationOperators; i++) {
			cluster_evo_data_x[group_id]._currentRelativeProb[seed][i] = cluster_evo_data_x[group_id]._parents[seed][cluster_evo_data_x[group_id]._currentParent[seed]][i];
		}
		// double * curr_parent = cluster_evo_data_x[group_id]._parents[seed][cluster_evo_data_x[group_id]._currentParent[seed]];
		// select a single random gene and invert
		for(int i =1; i <= 6; i++)
		{
			int randomGene = rand_32_int(kAflNumMutationOperators);
#ifdef REAL_VALUED
			double mutation = (rand_32_double_gauss() - 0.5) / 4 *20;
			cluster_evo_data_x[group_id]._currentRelativeProb[seed][randomGene] += mutation;
			if(cluster_evo_data_x[group_id]._currentRelativeProb[seed][randomGene] < 0)
				cluster_evo_data_x[group_id]._currentRelativeProb[seed][randomGene] = 0;
		}
		double *New_child = cluster_evo_data_x[group_id]._currentRelativeProb[seed];
		
		/*int randomGene = rand_32_int(kAflNumMutationOperators);
#ifdef REAL_VALUED
		double mutation = (rand_32_double_gauss() - 0.5) / 4;
		cluster_evo_data_x[group_id]._currentRelativeProb[seed][randomGene] += mutation;
		if(cluster_evo_data_x[group_id]._currentRelativeProb[seed][randomGene] < 0)
			cluster_evo_data_x[group_id]._currentRelativeProb[seed][randomGene] = 0;*/
#else		
		cluster_evo_data_x[group_id]._currentRelativeProb[seed][randomGene] = !cluster_evo_data_x[group_id]._currentRelativeProb[seed][randomGene];
#endif
	}
}

/*
 * Get the best parent solution so far as a vector (hardcoded to max mutation operators in AFL)
 * @seed: seed to attribute this to
*/
u32 MesFuzz_get_parent_repr(u32 seed, u32 group_id, cluster_data_structures_for_AMBES cluster_evo_data_x[]) {
	group_id -= 1;
	u32 value;
	bool *_currentParentBool = (bool *) (cluster_evo_data_x[group_id]._parents[seed][cluster_evo_data_x[group_id]._bestSoFar[seed]]); 
	for (int i = 0; i < 15; i++) {
		printf("%d: %s\n", i,_currentParentBool[i] ? "true" : "false");
	}
	/* encode parent in a number */
	for (int i = 14; i >= 0; i--) {
		value |= (cluster_evo_data_x[group_id]._parents[seed][cluster_evo_data_x[group_id]._bestSoFar[seed]][i]) ? (1 << i) : 0;
	}
	printf("val: %d\n", value);
	return value;
}



#define USE_FAST_RNG

#ifdef USE_FAST_RNG

#define ROTL(d,lrot) ((d<<(lrot)) | (d>>(8*sizeof(d)-(lrot))))
#define TWO31 0x80000000u
#define TWO32f (TWO31*2.0)

#endif


#define sigma 0.03
#define mu 0.5
#define two_pi 2*M_PI

bool initialized = false;

uint64_t xState = 0x0DDB1A5E5BAD5EEDull,
     yState = 0x519fb20ce6a199bbull;  // set to nonzero seed

void rand_init() {
    #ifdef USE_FAST_RNG
	if (!initialized){
	    srand(time(NULL));
            xState = (((uint64_t)(unsigned int)rand() << 32) + (uint64_t)(unsigned int)rand());
            yState = (((uint64_t)(unsigned int)rand() << 32) + (uint64_t)(unsigned int)rand());
	}
    #else
	if (!initialized) srand(time(NULL)):
    #endif
}

uint64_t romuDuoJr_random () {
    uint64_t xp = xState;
    xState = 15241094284759029579u * yState;
    yState = yState - xp;  yState = ROTL(yState,27);

    return xp;
}

inline double rand_32_double() {
    #ifdef USE_FAST_RNG
        unsigned int in = (unsigned int) romuDuoJr_random();
        double y = (double) in;
        return y/TWO32f;
    #else
        return (double) rand() / RAND_MAX;
    #endif
}

unsigned int rand_32_int (unsigned int max) {
    #ifdef USE_FAST_RNG
        unsigned int in = (unsigned int) romuDuoJr_random();
    #else
        unsigned int in = rand()
    #endif

    return in % max;
}

inline double rand_32_double_gauss() {
	double epsilon = DBL_EPSILON;

	double u1, u2;
	do
	{
		u1 = rand_32_double();
		u2 = rand_32_double();
	}
	while (u1 <= epsilon);

	double z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	// auto z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2); // not used here!
	//
	return fmin(fmax(z0 * sigma + mu, 0.0), 1.0);
}

inline double generate_uniform_random() {
    
    double sqrt12sigma = sqrt(12 * sigma);
    double a = mu - sqrt12sigma / 2;
    double b = mu + sqrt12sigma / 2;

    
    double u = (double)rand() / RAND_MAX;

    
    return a + u * (b - a);
}
