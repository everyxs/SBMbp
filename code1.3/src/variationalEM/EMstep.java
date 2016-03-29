package variationalEM;

import graphTools.Graph;
import java.util.Random;


/**
 * This class implements the abstract EMstep class
 * To be instantiated by BP and MCMC classes
 * Built for undirected multi-graphs and Poisson/DC block models
 * 
 * @author Xiaoran Yan ( everyxt@gmail.com )
 * @version BP_1.3
 * @time Nov, 2013
 */

public abstract class EMstep {
	
	// --- Instance Variables ----------------------------------------------------
	public boolean degreeCorrect; //flag for degree correction in model
	public boolean gSizeCorrect; //flag for group size correction in model
	public boolean fixEstep;
	public Graph graph;
	public int steps;
	double epsilon; //the threshold for convergence test
	public double[] gNode; //the group distribution of vertices
	public double[][] typeP; //the P_ij group affinity matrix
	static double[] factTable; //lookup table for factorial calculation

	public double likelihood;
	public double likelihoodHard;
	public double like3;
	public double theory;
	// --- Constructors ----------------------------------------------------------
	public EMstep(){} //the empty constructor
	/**
	 * This constructor creates a EMstep with random block parameters (to be overwritten)
	 * @param g Graph
	 * @param degreeC boolean
	 * @param gSize boolean
	 */
	public EMstep(Graph g, boolean degreeC, boolean gSize) {
		graph = g;
		degreeCorrect = degreeC;
		gSizeCorrect = gSize;
		fixEstep = false;
		steps = 0;
		epsilon = 0.0001 * graph.getNumNodes() * graph.getNumNodes() * graph.getNumType();
		gNode = new double[graph.getNumType()];
		typeP = new double[graph.getNumType()][graph.getNumType()]; 
		factTable = new double[g.getOutDegreeTop()+1]; //pending: memory bound
		factTable[0] = 1;
		for (int i=1; i<factTable.length; i++)
			factTable[i] = factTable[i-1] * i;
		likelihood = -Double.MAX_VALUE;
	}
	
	/**
	 * This constructor creates a EMstep with manually set block parameters, needs to be instantiated
	 * @param g Graph
	 * @param degreeC boolean
	 * @param gSize boolean
	 * @param p double[][]
	 * @param n double[]
	 */
	public EMstep(Graph g, boolean degree, boolean gSize, double[][] p, double[] n) {
		Random r = new Random();
		graph = g;
		degreeCorrect = degree;
		gSizeCorrect = gSize;
		fixEstep = false;
		steps = 0;
		epsilon = 0.0001 * graph.getNumNodes() * graph.getNumNodes() * graph.getNumType();
		gNode = new double[graph.getNumType()];
		for (int i=0; i<gNode.length; i++)
			gNode[i] = 1e-64;
		typeP = new double[graph.getNumType()][graph.getNumType()]; 
		for (int i=0; i<typeP.length; i++)
			for (int j=0; j<typeP[i].length; j++) {
				typeP[i][j] = r.nextDouble();
				if (degreeCorrect)
					typeP[i][j] = typeP[i][j] / java.lang.Math.pow(((double) graph.getNumEdgs() / graph.getNumNodes()), 2);
			}
		factTable = new double[g.getOutDegreeTop()+1]; //pending: memory bound
		factTable[0] = 1;
		for (int i=1; i<factTable.length; i++)
			factTable[i] = factTable[i-1] * i;
		likelihood = -Double.MAX_VALUE;
		update(p, n);
	}
	
	// --- Instance Methods ------------------------------------------------------
	/**
	 * This method calculates factorial. 
	 * @param t int
	 */	
	protected double poisson(double base, int pow) {
		if (pow==0)
			return 1;
		return Math.pow(base, pow) / factTable[pow];	
	}
	/**
	 * This method finds the maximum index in an array. 
	 * @param array double[]
	 */	
	protected int findMax(double[] array) {;
		int temp = 0;
		for (int i=1; i<array.length; i++)
			if (array[i]>array[temp])
				temp = i;
		return temp;	
	}
	/**
	 * This method updates the hyper parameters after the maximization step.
	 * @param p double[][]
	 * @param n double[]
	 */
	public void update(double[][] p, double[] n) {
		steps = 0;
		likelihood = 0;
		for (int i=0; i<p.length; i++) {
			for (int j=0; j<p[0].length; j++)
				typeP[i][j] = p[i][j];
			gNode[i] = n[i];
		}
	}
	/**
	 * This method updates the hyper parameters (abstract method to be instantiated).
	 */
	public abstract double mStep(double[] n, boolean fix);
	
	/**
	 * This method gives a random permutation of nodes for message[][].
	 * @return a random list for message updates
	 * @param null
	 */
	public int[] permute () {
		Random r = new Random();
		int[] updateList = new int[graph.getNumNodes()];
		for (int i=0; i<updateList.length; i++)
			updateList[i] = i;
		for (int i=0; i<updateList.length; i++) {
			int j = r.nextInt(updateList.length);
			int temp = updateList[i];
			updateList[i] = updateList[j];
			updateList[j] = temp;
		}
		return updateList;
	}

	/**
	 * This method repeats the flooding steps until convergence,
	 * and returns the stationary marginals (abstract method to be instantiated).
	 */
	public abstract double[][] convergeExpectation();
	
}