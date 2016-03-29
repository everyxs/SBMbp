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

public class MCMC extends EMstep{
	
	// --- Instance Variables ----------------------------------------------------
	public int count; //count the number of samples across initials
	double beta; //the parameter correspond to the inverse temperature in statistics physics
	double[][] marginal; //the accumulated marginal mixture for each vertex
	double[][] convergeTest; //marginal mixture for convergence testing in each initial
	double[][] tP; //temporary p_ij
	double meanLogLike;
	double maxLogLike;
	private Classification classification; //subclass for MCMC sampling
	// --- Constructors ---------------------------------------------------------- 
	public MCMC(){} //the empty constructor
	/**
	 * This constructor creates a MCMC sampler with manually set block parameters and specified temperature
	 * @param g Graph
	 * @param degreeC boolean
	 * @param gSize boolean
	 * @param typeP double[][]
	 * @param gNode double[]
	 * @param invTemp double
	 */
	public MCMC(Graph g, boolean degreeC, boolean gSize, double[][] typeP, double[] gNode, double invTemp) {
		super(g, degreeC, gSize, typeP, gNode);
		epsilon = 0.00001 * graph.getNumNodes() * graph.getNumType();
		count = graph.getNumNodes()*graph.getNumType()*5;
		beta = invTemp;
		marginal = new double[graph.getNumNodes()][graph.getNumType()];
		convergeTest = new double[graph.getNumNodes()][graph.getNumType()];
		tP = new double[graph.getNumType()][graph.getNumType()];
		meanLogLike = 0;
		maxLogLike = -Double.MAX_VALUE;
		classification = new Classification(graph, graph.getNumType());
	}
	/**
	 * This constructor creates a MCMC sampler from an exact copy, keeping the initial classifications
	 * @param g Graph
	 * @param copy MCMC
	 * @param invTemp double
	 */
	public MCMC(Graph g, MCMC copy, double invTemp) {
		super(g, copy.degreeCorrect, copy.gSizeCorrect, copy.typeP, copy.gNode);
		tP = new double[graph.getNumType()][graph.getNumType()];
		convergeTest = new double[graph.getNumNodes()][graph.getNumType()];
		beta = invTemp;
		
		epsilon = copy.epsilon;
		count = copy.count;
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
			distribution = classification.distribution(Ulist[i], beta, gSizeCorrect, degreeCorrect, typeP, gNode); //get the heat bath MCMC distribution
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
		maxLogLike = -Double.MAX_VALUE;
		// update in order according to the permutation
		Random r = new Random();
		double[] distribution;
		
		for (int i=0; i<Ulist.length; i++) {
			distribution = classification.distribution(Ulist[i], beta, gSizeCorrect, degreeCorrect, typeP, gNode); //get the heat bath MCMC distribution
			for (int j = 0; j < classification.numGroup(); j++) {
				marginal[Ulist[i]][j] +=  distribution[j];
				for (int k=0; k<classification.numGroup(); k++)
					tP[j][k] += classification.avgGroupMatrix[0][j][k];
			}
			//Change the classification according to the distribution
			double randgroup = r.nextDouble();//[0d,1d)
			int group = -1;
			while (randgroup >= 0) {
				group++;
				randgroup = randgroup - distribution[group];
			}
			double temp = classification.likelihood(1, gSizeCorrect, degreeCorrect, typeP, gNode);
			meanLogLike += temp;
			if (temp > maxLogLike)
				maxLogLike = temp;
			classification.mutate(Ulist[i], group);
		}
	}
	
	/**
	  * This method implements the M-step after the E-step converges
	  * A polymorphic extension of the abstract method mStep in the parent class
	  * Built for undirected multi-graphs and Poisson/DC block models
	  * @return maximum likelihood after this EM iteration
	  * @param marginals double[]
	  * @param fix boolean
	 */
	public double mStep(double[] marginals, boolean fix) {
		
		double delta = 0; //change measure
		for (int k1=0; k1<graph.getNumType(); k1++) {
			if (gSizeCorrect) //group size correction
				for (int i=0; i<graph.getNumNodes(); i++)
					likelihood += java.lang.Math.log(marginal[i][k1]); // left shifted to avoid overflow
			gNode[k1] = marginals[k1];
			for (int k2=0; k2<graph.getNumType(); k2++) {
				tP[k1][k2] = tP[k1][k2] / graph.getNumNodes() / count;
				delta += java.lang.Math.abs(tP[k1][k2]-typeP[k1][k2]);
				typeP[k1][k2] = tP[k1][k2];
			}
		}
			
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
		
		marginal = new double[graph.getNumNodes()][graph.getNumType()];
		tP = new double[graph.getNumType()][graph.getNumType()];
		meanLogLike = 0;
		
		convergeTest = new double[graph.getNumNodes()][graph.getNumType()];
		int[] topList = new int[1];
		topList[0] = -1;
		steps = 0;
		double delta = Double.MAX_VALUE; //iteration level difference
		//the burn-in steps until convergence
		while (delta > epsilon && steps<5) {
			delta = stepMC();
			steps++;
		}
		//statistics collecting after the MC is well mixed
		for (int i=0; i<count; i++)
			stepMCmixed();
		
		//likelihood  = meanLogLike / graph.getNumNodes() / count;
		likelihood = maxLogLike;
		for (int i=0; i<marginal.length; i++)
			for (int j=0; j<marginal[0].length; j++)
				marginal[i][j] = marginal[i][j] / count;
		
		return marginal;
	}
	
}