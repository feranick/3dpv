/***************************************************************/
 /* Single & Multi-Objective Real-Coded Genetic Algorithms Code */
/* Author: Kumara Sastry                                       */
/* Illinois Genetic Algorithms Laboratory (IlliGAL)            */
/* Deparment of General Engineering                            */
/* University of Illinois at Urbana-Champaign                  */
/* 104 S. Mathews Ave, Urbana, IL 61801                        */
/***************************************************************/

#include "crossover.hpp"
// All constructors first - pretty straightforward stuff
OneTwoPointCrossover::OneTwoPointCrossover(int numPoints) :
	noOfPoints(numPoints) {
}

UniformCrossover::UniformCrossover(void) :
	genewiseProbability(0.5) {
}

UniformCrossover::UniformCrossover(double geneProb) :
	genewiseProbability(geneProb) {
}

SimulatedBinaryCrossover::SimulatedBinaryCrossover(void) :
	genewiseProbability(0.5) {
}
SimulatedBinaryCrossover::SimulatedBinaryCrossover(double geneProb) :
	genewiseProbability(geneProb) {
}

/** 
 Single and Two Point Crossover: 
 Single Point: Randomly select a crossover point and exchange the genes to the left of the crossover point to create two new individuals.
 Two Point: Randomly select two crossover points and exchange genes between the two crossover points to creat two new indivduals.
 */
void OneTwoPointCrossover::crossover(Individual *parent1, Individual *parent2) {
	int ii, XoverPt1, XoverPt2;
	double Temp;
	Individual child1(parent1), child2(parent2);

	switch (noOfPoints) {
	case 1:
		///////////DEFAULT ONEPOINT CROSS
		/*XoverPt1 = myRandom.boundedIntegerRandom(0,
				globalSetup->noOfDecisionVariables-1);
		for (ii = 0; ii <= XoverPt1; ii++) {
			Temp = child1[ii];
			child1.setValue(ii, child2[ii]);
			child2.setValue(ii, Temp);
		}*/
		//////////////////TRIANGLE SAWP METHOD FOR 9 decision VARS per CELL////////////////
		XoverPt1 = myRandom.boundedIntegerRandom(0, globalSetup->noOfDecisionVariables-1);
		XoverPt1 = XoverPt1 - (XoverPt1%9);
		// originally the parents only switch genes to the right of this point, but that seems limiting
		// so I will make a 1/2 probability of either swapping the left chunk or right chunk
		if (myRandom.boundedIntegerRandom(0,1))
		{
			// swap genes to left of crosspoint
			for (ii = 0; ii <= XoverPt1; ii++) {
				Temp = child1[ii];
				child1.setValue(ii, child2[ii]);
				child2.setValue(ii, Temp);
			}
		}
		else
		{
			// swap genes to the right of crosspoint
			for (ii = XoverPt1; ii < globalSetup->noOfDecisionVariables; ii++) {
				Temp = child1[ii];
				child1.setValue(ii, child2[ii]);
				child2.setValue(ii, Temp);
			}
		}
		break;
	case 2:
		//////////DEFAULT//////////////////////////
		/*XoverPt1 = myRandom.boundedIntegerRandom(0,
				globalSetup->noOfDecisionVariables-1);
		do {
			XoverPt2 = myRandom.boundedIntegerRandom(0,
					globalSetup->noOfDecisionVariables-1);
		} while (XoverPt1 == XoverPt2);
		if (XoverPt1 > XoverPt2)
			SWAP(XoverPt1, XoverPt2);
		for (ii = XoverPt1; ii < XoverPt2; ii++) {
			Temp = child1[ii];
			child1.setValue(ii, child2[ii]);
			child2.setValue(ii, Temp);
		}*/
		/////TRIANGLE SWAP METHOD FOR 9 decision VARS per cell///////////////////////
		XoverPt1 = myRandom.boundedIntegerRandom(0, globalSetup->noOfDecisionVariables-1);
		XoverPt1 = XoverPt1 - (XoverPt1%9);   // get a point = 0mod9 which is always start of a triangle list of 9 coordinates
		do {
			XoverPt2 = myRandom.boundedIntegerRandom(0, globalSetup->noOfDecisionVariables-1);
			XoverPt2 = XoverPt2 - (XoverPt2%9);   // get a point = 0mod9 which is always start of a triangle list of 9 coordinates
		} while (XoverPt1 == XoverPt2);
		if (XoverPt1 > XoverPt2)
			SWAP(XoverPt1, XoverPt2);
		XoverPt2+=1; // w/o this it seems the last coordinate of the last triangle would be left out of cross
					 // so I'm defining point 2 to be included in the crossing
		for (ii = XoverPt1; ii < XoverPt2; ii++) {   // < XoverPt2 ensures last switched point is a last coordinate in one triangle
			Temp = child1[ii];
			child1.setValue(ii, child2[ii]);
			child2.setValue(ii, Temp);
		}
		break;
	default:
		exit(0);
	}
	if (globalSetup->nichingType == DeterministicCrowding) {
		child1.evaluateFitness();
		child2.evaluateFitness();
		deterministicCrowding(parent1, parent2, &child1, &child2);
	} else {
		for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++) {
			parent1->setValue(ii, child1[ii]);
			parent2->setValue(ii, child2[ii]);
		}
	}
}

/// Uniform crossover: swap each gene with a probability of genewiseProbability (default value of 0.5)
// BAM: I do not think this would be good to use in the triangle-coordinate picture of the decision variables
// BAM: (3/13/09) - maybe ill try this again.....now that my program is more mature

void UniformCrossover::crossover(Individual *parent1, Individual *parent2) {
	int ii;
	double Temp;
	Individual child1(parent1), child2(parent2);

	if (myRandom.flip(globalSetup->xOverProbability)) {
		for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++) {
			if (myRandom.flip(genewiseProbability)) {
				Temp = child1[ii];
				child1.setValue(ii, child2[ii]);
				child2.setValue(ii, Temp);
			}
		}
		if (globalSetup->nichingType == DeterministicCrowding) {
			child1.evaluateFitness();
			child2.evaluateFitness();
			deterministicCrowding(parent1, parent2, &child1, &child2);
		} else {
			for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++) {
				parent1->setValue(ii, child1[ii]);
				parent2->setValue(ii, child2[ii]);
			}
		}
	}
}

/**
 Simulated Binary Crossover: Creates new individuals based on a probability distribution (a polynomial one) similar to binary coded single point crossover. 
 Reference: Deb, K., & Agarwal, R.B. (1995) "Simulated Binary Crossover for Continuous Search Space", Complex Systems. 9. pp. 115-148. (TCGA No. 05896).
 Deb, K., & Kumar, A. (1995) "Real-Coded Genetic Algorithms with Simulated Binary Crossover: Studies on Multimodal and Multiobjective Problems", Complex Systems. 9. pp. 431-454. (TCGA No. 05897).
 */
void SimulatedBinaryCrossover::crossover(Individual *parent1,
		Individual *parent2) {
	int ii;
	double gene1Value, gene2Value, beta, betaq, geneMin, geneMax, alpha;
	double tempGeneValue;
	double etaC; //The order of the polynomial for the SBX crossover
	Individual child1(parent1), child2(parent2);

	if (globalSetup->xOverParameters==NULL)
		etaC = 10;
	else
		etaC = ((double *)globalSetup->xOverParameters)[1];
	if (myRandom.flip(globalSetup->xOverProbability)) {
		for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++) {
			if ((myRandom.flip(genewiseProbability)) && (fabs(child1[ii]
					-child2[ii]) >= ZERO)) {
				geneMin = globalSetup->variableRanges[ii][0];
				geneMax = globalSetup->variableRanges[ii][1];
				if (child2[ii] > child1[ii]) {
					gene1Value = child1[ii];
					gene2Value = child2[ii];
				} else {
					gene1Value = child2[ii];
					gene2Value = child1[ii];
				}
				if ((gene1Value - geneMin) > (geneMax - gene2Value))
					beta = (gene2Value-gene1Value)/(2.0*geneMax - gene2Value
							- gene1Value);
				else
					beta = (gene2Value-gene1Value)/(gene2Value + gene1Value
							- 2.0*geneMin);
				alpha = myRandom.random01()*(2.0 - pow(beta, (double)(etaC+1)));
				if (alpha <= 1.0)
					betaq = pow(alpha, 1.0/((double)(etaC+1)));
				else
					betaq = pow(1.0/(2.0-alpha), 1.0/((double)(etaC+1)));

				tempGeneValue = 0.5*((gene1Value+gene2Value) - betaq
						*(gene2Value-gene1Value));
				if (globalSetup->variableTypes[ii] == Integer)
					tempGeneValue = int(tempGeneValue+0.5);
				child1.setValue(ii, tempGeneValue);
				tempGeneValue = 0.5*((gene1Value+gene2Value) + betaq
						*(gene2Value-gene1Value));
				if (globalSetup->variableTypes[ii] == Integer)
					tempGeneValue = int(tempGeneValue+0.5);
				child2.setValue(ii, tempGeneValue);
			}
		}
		if (globalSetup->nichingType == DeterministicCrowding) {
			child1.evaluateFitness();
			child2.evaluateFitness();
			deterministicCrowding(parent1, parent2, &child1, &child2);
		} else {
			for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++) {
				parent1->setValue(ii, child1[ii]);
				parent2->setValue(ii, child2[ii]);
			}
		}
	}
}

/**
 Deterministic Crowding: Compare each children obtained through crossover with the nearest parent based on phenotypic distance and replace the child with the nearest parent if the parent has higher fitness
 Reference: Mahfoud, S.W. (1992) "Crowding and Preselection Revisited", In R. Manner and B. Manderic (Eds.) Parallel Problem Solving From Nature, 2. (pp. 27-36). (IlliGAL Report No. 92004). http://www-illigal.ge.uiuc.edu/techreps.php3. (ftp://ftp-illigal.ge.uiuc.edu/pub/papers/IlliGALs/92004.ps.Z).
 */
void Crossover::deterministicCrowding(Individual *parent1, Individual *parent2,
		Individual *child1, Individual *child2) {

	double phenotypeDist[4] = { 0.0, 0.0, 0.0, 0.0 };
	int ii;

	for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++) {
		phenotypeDist[0] += (((*parent1)[ii]-(*child1)[ii])*((*parent1)[ii]-(*child1)[ii]));
		phenotypeDist[1] += (((*parent1)[ii]-(*child2)[ii])*((*parent1)[ii]-(*child2)[ii]));
		phenotypeDist[2] += (((*parent2)[ii]-(*child1)[ii])*((*parent2)[ii]-(*child1)[ii]));
		phenotypeDist[3] += (((*parent2)[ii]-(*child2)[ii])*((*parent2)[ii]-(*child2)[ii]));
	}

	if (phenotypeDist[0]+phenotypeDist[3] <= phenotypeDist[1]+phenotypeDist[2]) {
		if (isBetter(child1, parent1))
			for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++)
				parent1->setValue(ii, (*child1)[ii]);
		if (isBetter(child2, parent2))
			for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++)
				parent2->setValue(ii, (*child2)[ii]);
	} else {
		if (isBetter(child2, parent1))
			for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++)
				parent1->setValue(ii, (*child2)[ii]);
		if (isBetter(child1, parent2))
			for (ii = 0; ii < globalSetup->noOfDecisionVariables; ii++)
				parent2->setValue(ii, (*child1)[ii]);
	}
}

