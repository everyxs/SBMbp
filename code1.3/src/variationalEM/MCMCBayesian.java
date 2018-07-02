package variationalEM;

import graphTools.Graph;

import java.util.Random;

/**
 * This class implements the MAP MCMC for comparison
 * A child class extends the abstract EMstep class
 * Built for undirected multi-graphs and Poisson/DC block models
 * 
 * @author Xiaoran Yan ( everyxt@gmail.com )
 * @version BP_1.3
 * @time Nov, 2013
 */

public class MCMCBayesian extends EMstep{
	
	// --- Instance Variables ----------------------------------------------------
	public int count; //count the number of samples across initials
	int groupCount; //number of blocks in SBM
	double beta; //the parameter correspond to the inverse temperature in statistics physics
	double[][] marginal; //the accumulated marginal mixture for each vertex
	double[][] convergeTest; //marginal mixture for convergence testing in each initial
	double[][] tP; //temporary p_ij
	double meanLogLike;
	double maxLogLike;
	private ClassificationBayesian classification; //subclass for MCMC sampling
	// --- Constructors ---------------------------------------------------------- 
	public MCMCBayesian(){} //the empty constructor
	/**
	 * This constructor creates a MCMC sampler with manually set block parameters and specified temperature
	 * @param g Graph
	 * @param degreeC boolean
	 * @param gSize boolean
	 * @param typeP double[][]
	 * @param gNode double[]
	 * @param invTemp double
	 */
	public MCMCBayesian(Graph g, boolean degreeC, double[][] typeP, double[] gNode, double invTemp, int k) {
		super(g, degreeC, true, typeP, gNode);
		epsilon = 0.0000001 * graph.getNumNodes() * graph.getNumType();
		count = graph.getNumNodes()*k*2;
		groupCount = k;
		beta = invTemp;
		marginal = new double[graph.getNumNodes()][k];
		convergeTest = new double[graph.getNumNodes()][k];
		tP = new double[k][k];
		meanLogLike = 0;
		maxLogLike = -Double.MAX_VALUE;
		classification = new ClassificationBayesian(graph, groupCount, true);
	}
	/**
	 * This constructor creates a MCMC sampler from an exact copy, keeping the initial classifications
	 * @param g Graph
	 * @param copy MCMC
	 * @param invTemp double
	 */
	public MCMCBayesian(Graph g, MCMCBayesian copy, double invTemp) {
		super(g, copy.degreeCorrect, copy.gSizeCorrect, copy.typeP, copy.gNode);
		tP = new double[copy.groupCount][copy.groupCount];
		convergeTest = new double[graph.getNumNodes()][copy.groupCount];
		beta = invTemp;
		
		epsilon = copy.epsilon;
		count = copy.count;
		groupCount = copy.groupCount;
		marginal = copy.marginal;
		likelihood = copy.likelihood;
		meanLogLike = copy.meanLogLike;
		maxLogLike = copy.maxLogLike;
		classification = copy.classification;
		
	}
	// --- Instance Methods ------------------------------------------------------
	/**
	 * This method does a single step of MCMC update (asynchronous) across the network,
	 * This is the converging version which returns a measure of change in l1 norm
	 * @param null
	 */
	public double stepMC() {
		int[] Ulist = permute();
		double delta = 0;
		
		// update in order according to the permutation
		Random r = new Random();
		double[] distribution;
		
		for (int i=0; i<Ulist.length; i++) {
			distribution = classification.distribution(Ulist[i], beta, gSizeCorrect, degreeCorrect); //get the heat bath MCMC distribution
			delta = 0;
			for (int j = 0; j < classification.numGroup(); j++) {
				convergeTest[Ulist[i]][j] +=  distribution[j];
				delta += java.lang.Math.abs(distribution[j] - convergeTest[Ulist[i]][j] / steps);
			}
			//Change the classification according to the distribution
			double randgroup = r.nextDouble();//[0d,1d)
			int group = -1;
			while (randgroup >= 0) {
				group++;
				randgroup = randgroup - distribution[group];
			}
			classification.mutate(Ulist[i], group);
		}
		return delta;
	}
	
	/**
	 * This method does a single step of MCMC update (asynchronous) across the network,
	 * This is the mixed version for statistics collection
	 * @param null
	 */
	public void stepMCmixed() {
		int[] Ulist = permute();
		// update in order according to the permutation
		Random r = new Random();
		double[] distribution;
		
		for (int i=0; i<Ulist.length; i++) {
			distribution = classification.distribution(Ulist[i], beta, gSizeCorrect, degreeCorrect); //get the heat bath MCMC distribution
			for (int j = 0; j < classification.numGroup(); j++) {
				marginal[Ulist[i]][j] +=  distribution[j];
			}
			//Change the classification according to the distribution
			double randgroup = r.nextDouble();//[0d,1d)
			int group = -1;
			while (randgroup >= 0) {
				group++;
				randgroup = randgroup - distribution[group];
			}
			double temp = classification.likelihood(beta, gSizeCorrect, degreeCorrect);
			meanLogLike += temp;
			if (temp > maxLogLike)
				maxLogLike = temp;
			classification.mutate(Ulist[i], group);
		}
	}
	
	/**
	  * This method implements(dummy) the M-step after the E-step converges
	  * For MAP optimization, the M-step is skipped
	  * @return maximum likelihood after this EM iteration
	  * @param marginals double[]
	  * @param fix boolean
	 */
	public double mStep(double[] marginals, boolean fix) {
		
		double delta = 0; //change measure
		return delta;
	}
	
	/**
	  * This method implements the E-step inner-loops for the MCMC sampler
	  * A polymorphic extension of the abstract method convergeExpectation in parent class
	  * Built for undirected multi-graphs and Poisson/DC block models
	  * @return block marginal vectors for all nodes
	  * @param null
	 */
	public double[][] convergeExpectation() {
		
		marginal = new double[graph.getNumNodes()][groupCount];
		tP = new double[graph.getNumType()][groupCount];
		meanLogLike = 0;
		//maxLogLike = -Double.MAX_VALUE;
		
		convergeTest = new double[graph.getNumNodes()][groupCount];
		int[] topList = new int[1];
		topList[0] = -1;
		steps = 0;
		
		//the burn-in steps until convergence
		//while (delta > epsilon && steps<count) {
		//	delta = stepMC();
		//	steps++;
		//}
		//statistics collecting after the MC is well mixed
		//if (steps < count/10)
			//steps = count/10;
		for (int i=0; i<count; i++)
			stepMCmixed();
		
		likelihood = maxLogLike;
		marginal[0][0] = steps;
		marginal[1][0] = count;
		classification = new ClassificationBayesian(classification); //reinitialize MCMC for next outer iteration
		return marginal; //dummy output
	}
	
}