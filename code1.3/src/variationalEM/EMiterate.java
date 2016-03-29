package variationalEM;

import java.util.Random;

/**
 * This class implements the EM outer-loop using EMstep as the inner-loops
 * 
 * 
 * @author Xiaoran Yan ( everyxt@gmail.com )
 * @version BP_1.3
 * @time Nov, 2013
 */

public class EMiterate {
	// --- Instance Variables ----------------------------------------------------
	EMstep emStep; //the abstract class of emStep
	int iterations; //number of iterations in the outer-loop
	double epsilon; //convergence threshold
	public double[][] typeP; //the P_ij group affinity matrix
	public double[][] marginal; //the mixed membership vector
	
	// --- Constructors ----------------------------------------------------------
	public EMiterate(){} //the empty constructor
	/**
	 * This constructor initializes an EM process.
	 * @param em EMstep
	 * @param assort double 
	 * @param disassort double
	 */
	public EMiterate(EMstep em, double assort, double disassort) {
		iterations = 0;
		emStep = em;
		//The P_ij group affinity matrix
		typeP = new double[emStep.graph.getNumType()][emStep.graph.getNumType()]; 
		for (int i=0; i<typeP.length; i++)
			for (int j=0; j<typeP[0].length; j++){
				if (i==j)
					typeP[i][j] = assort;
				else
					typeP[i][j] = disassort;
			}
		
		epsilon = 0.001 * emStep.graph.getNumEdgs() / java.lang.Math.pow(emStep.graph.getNumNodes(),2);
		
		
		// Correct P_ij by degree
		if (emStep.degreeCorrect) {
			for (int i=0; i<typeP.length; i++)
				for (int j=0; j<typeP[0].length; j++){
					typeP[i][j] = typeP[i][j]/4 / java.lang.Math.pow(((double) emStep.graph.getNumEdgs() / emStep.graph.getNumNodes()), 2);
				}
			epsilon = epsilon / java.lang.Math.pow(((double) emStep.graph.getNumEdgs() / emStep.graph.getNumNodes()), 2);
		}
		
		
		//The mixed membership vector
		marginal = new double[emStep.graph.getNumNodes()][emStep.graph.getNumType()];
		for (int i=0; i<marginal.length; i++)
			for (int j=0; j<marginal[0].length; j++){
				marginal[i][j] = 1.0 / emStep.graph.getNumType(); //uniform initialization
			}
	}
	
	/**
	 * This constructor initializes an EM process to manually set parameters (for parameter inheritance in tuning)
	 * @param em EMstep
	 * @param p double[][]
	 * @param M double[][]
	 * @param tune boolean
	 */
	public EMiterate(EMstep em, double[][] P, double[][] M, boolean tune) {
		iterations = 0;
		emStep = em;
		//The P_ij group affinity matrix
		typeP = new double[emStep.graph.getNumType()][emStep.graph.getNumType()]; 
		for (int i=0; i<P.length; i++)
			for (int j=0; j<P[0].length; j++)
				typeP[i][j] = P[i][j];
				
		epsilon = 0.001 * emStep.graph.getNumType() * emStep.graph.getNumType();
		
		if (emStep.degreeCorrect) //mean-field degree scaling
			epsilon = epsilon / java.lang.Math.pow(((double) emStep.graph.getNumEdgs() / emStep.graph.getNumNodes()), 2);
		
		//The mixed membership vector
		marginal = new double[emStep.graph.getNumNodes()][emStep.graph.getNumType()];
		for (int i=0; i<M.length; i++)
			for (int j=0; j<M[i].length; j++){
				marginal[i][j] = M[i][j];
			}
	}
	/**
	 * This constructor initializes an EM process from a copy
	 * @param em EMstep
	 * @param copy EMiterate
	 */
	public EMiterate(EMstep em, EMiterate copy) {
		iterations = 0;
		emStep = em;
		//The P_ij group affinity matrix
		typeP = new double[emStep.graph.getNumType()][emStep.graph.getNumType()]; 
		for (int i=0; i<typeP.length; i++)
			for (int j=0; j<typeP[0].length; j++)
				typeP[i][j] = copy.typeP[i][j];
				
		epsilon = copy.epsilon;
		//The mixed membership vector
		marginal = new double[emStep.graph.getNumNodes()][emStep.graph.getNumType()];
		for (int i=0; i<marginal.length; i++)
			for (int j=0; j<marginal[0].length; j++){
				marginal[i][j] = copy.marginal[i][j];
			}
	}
	
	// --- Instance Methods ------------------------------------------------------
	/**
	 * This method calls a single iteration of EM, with the option to fix the M step.
	 * @param fix boolean
	 */
	public double stepEM(boolean fix) {
		//Expectation step
		marginal = emStep.convergeExpectation();
		
		//Maximization step
		double[] gNode = new double[emStep.graph.getNumType()];
		double sum= 0;
		for (int i=0; i<gNode.length; i++){
			for (int j=0; j<marginal.length;j++)
				gNode[i] += marginal[j][i];
			sum += gNode[i];
		}
		for (int i=0; i<gNode.length; i++)
			gNode[i] = gNode[i] / sum; //MLE for the node parameters
		double delta = emStep.mStep(gNode, fix); //switch for fixing M step
		for (int i=0; i<typeP.length; i++)
			for (int j=0; j<typeP.length; j++)
				typeP[i][j] = emStep.typeP[i][j]; //MLE for the edge parameters
		return delta; //return the change in parameter values 
	}
	
	/**
	 * This method repeats the EM steps until convergence (optional manual initialization)
	 * @return the maximum likelihood value double
	 * @param gN double[]
	 * @param tP double[][]
	 */
	public double convergeEM(double[] gN, double[][] tP) {
		double delta = Double.MAX_VALUE;

		Random r = new Random(); //default random initialization
		double[][] typeP = new double[emStep.graph.getNumType()][emStep.graph.getNumType()];
		for (int i=0; i<typeP.length; i++)
			for (int j=0; j<typeP[i].length; j++) {
				typeP[i][j] = r.nextDouble();
				if (emStep.degreeCorrect) //for DC model
					typeP[i][j] = typeP[i][j] / Math.pow(((double) emStep.graph.getNumEdgs() / emStep.graph.getNumNodes()), 2);
			}
		for (int i=0; i<tP.length; i++)
			for (int j=0; j<tP[i].length; j++)
				typeP[i][j] = tP[i][j]; //overwrite with manual initializations

		double[] gNode = new double[emStep.graph.getNumType()];
		for (int i=0; i<gNode.length; i++)
			gNode[i] = 1e-6; //default epsilon initialization
		for (int i=0; i<gN.length; i++)
			gNode[i] = gN[i]; //overwrite with manual initializations
		
		double sum = 0;
		for (int i=0; i<gNode.length; i++)
			sum += gNode[i];
		
		for (int i=0; i<gNode.length; i++)
			gNode[i] = gNode[i] / sum; //normalization
		
		emStep.update(typeP, gNode); //update the parameter values
		iterations= 0;
		while(delta>epsilon && iterations<5) { //bound the total iteration to 5
			delta = stepEM(false); //iterate EM steps
			iterations++;
			System.out.println("EM tune iteration "+iterations+" BP steps used: "+emStep.steps); //print the number of loops
		}
		return emStep.likelihood;
	}
	
	/**
	 * This method calls a single iteration of stepEM, with specified parameter values
	 * @return the maximum likelihood value double
	 * @param gN double[]
	 * @param tP double[][]
	*/
	public double convergeEMfix(double[] gN, double[][] tP) { 
		double[][] typeP = new double[emStep.graph.getNumType()][emStep.graph.getNumType()];
		for (int i=0; i<typeP.length; i++)
			for (int j=0; j<typeP[i].length; j++) { //null model for the edge parameters
				typeP[i][j] = Math.pow(((double) emStep.graph.getNumEdgs() / emStep.graph.getNumNodes()), 2);
			}
		for (int i=0; i<tP.length; i++)
			for (int j=0; j<tP[i].length; j++)
				if (tP[i][j]>0)
					typeP[i][j] = tP[i][j]; //setting the edge parameters

		double[] gNode = new double[emStep.graph.getNumType()];
		for (int i=0; i<gNode.length; i++)
			gNode[i] = 1e-6;
		for (int i=0; i<gN.length; i++)
			if (gN[i]>0)
				gNode[i] = gN[i]; //setting the node parameters
		double sum = 0;
		for (int i=0; i<gNode.length; i++)
			sum += gNode[i];
		for (int i=0; i<gNode.length; i++)
			gNode[i] = gNode[i] / sum; //normalization
		
		emStep.update(typeP, gNode); //update the parameter values
		iterations= -1; //indicates fixed M-step
		stepEM(true); //iterate a single EM step with the M-step fixed
		System.out.println("EM fix iteration "+iterations+" BP steps used: "+emStep.steps); //print the number of loops
		return emStep.likelihood;	
	} 
	
	/**
	 * This method repeats the stepEM until convergence (used in grid search without parameter specification)
	 * @return the maximum likelihood value double
	 * @param null
	 */
	public double convergeEMgrid() {
		
		double delta = Double.MAX_VALUE;
		// Initial iteration
		double[] gNode = new double[emStep.graph.getNumType()];
		for (int i=0; i<gNode.length; i++)
			for (int j=0; j<marginal.length;j++)
				gNode[i] += marginal[j][i];	
		emStep.update(typeP, gNode); //update the parameter values
		
		while(delta>epsilon && iterations<5) { //bound the total iteration to 10
			delta = stepEM(false); //iterate EM steps
			iterations++;
			System.out.println("EM restart iteration "+iterations+" BP steps used: "+emStep.steps); //print the number of loops
		}
		return emStep.likelihood;
	}
}