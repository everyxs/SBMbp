package variationalEM;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.Random;
import graphTools.Graph;

/**
 * This class includes the parameters and initializations for using the variational EM framework 
 * for Bayesian model selection of number of blocks and between the vanilla block model and degree corrected (DC) block model.
 * This example showcases how different number of blocks can be fitted to the same data sets under this framework, and
 * how the polymorphism mechanism can be used to employ different EM implementations.
 * Please only edit the main function at the bottom for testing purposes. 
 * For other uses of the variational EM framework, users can write their own main functions according the API.
 * 
 * @author Xiaoran Yan ( everyxt@gmail.com )
 * @version BP_1.3
 * @time Nov, 2013
 */
public class MainBayesian {
	// --- Instance Variables ----------------------------------------------------
	public Graph graph; //input graph
	public EMstep MCMC; //the maximum likelihood EM setup (abstract class to be instantiated) 
	public double likelihood; //the maximum likelihood
	public double likelihoodVanilla; //the maximum likelihood for vanilla model
	
	// --- Constructors ---------------------------------------------------------- 
	public MainBayesian() {} //the empty constructor
	
	/**
	 * This constructor initializes internal parameters given the input
	 * @param dir String
	 */
	public MainBayesian(String dir) {
		likelihood = -Double.MAX_VALUE;
		likelihoodVanilla = -Double.MAX_VALUE;
		
		FileReader input = null;
		try {
			input = new FileReader(dir);
		}
		catch (FileNotFoundException ex) {
			ex.printStackTrace();
		}
		//Create a graph with the input file
		graph = new Graph(input, false);
	}
	
	// --- Instance Methods ------------------------------------------------------
	/**
	 * This method fits both the vanilla block model and degree corrected (DC) block model to a sample graph
	 * Initial parameters are grid searched with an option to start around the ground truth for testing purposes
	 * @return likelihood ratio between the models double
	 * @param print PrintStream
	 * @param em EMstep (abstract class to be instantiated)
	 */
	public double likeRatio(PrintStream print, EMstep em) {
		
		Graph newGraph = graph.RandomizeEdge(false, false); //generate a new random graph from the vanilla block model
		print.print(newGraph.getNumNodes() + "\t");
		print.print(newGraph.getNumEdgs() + "\t");
		
		print.print("vanilla" + "\t"); //initializing parameters around the ground truth parameters (unused for grid search)
		posterior(newGraph, em); //grid search parameters values for maximum likelihood
		likelihoodVanilla = likelihood; //remember the maximum likelihood
		for (int i=0; i<newGraph.getNumType(); i++) { //print the final parameters
			print.print("["+MCMC.gNode[i]+"]"+"\t");
			for (int j=0; j<newGraph.getNumType(); j++) {
				print.print(MCMC.typeP[i][j]+"\t");				
			}
		}
		print.print(likelihood+"\t");
		
		// Fitting the degree corrected (DC) block model
		em.degreeCorrect = true;
		print.print("degree corrected" + "\t");
		posterior(newGraph, em); //grid search parameters values for maximum likelihood
		for (int i=0; i<graph.getNumType(); i++) { //print the final parameters
			print.print("["+MCMC.gNode[i]+"]"+"\t");
			for (int j=0; j<graph.getNumType(); j++) {
				print.print(MCMC.typeP[i][j]+"\t");
			}
		}
		print.print(likelihood+"\t");
		return likelihood - likelihoodVanilla; //return the log-likelihood difference
	}
	
	/**
	 * This method does a grid search for parameters values based on the EM framework and 
	 * calls the convergeEMgrid method which requires no specified parameter inputs
	 * Use this method when there is no domain knowledge about the block structures
	 * @param g Graph
	 * @param em EMstep (abstract class to be instantiated)
	 */
	public double[][] posterior(Graph g, EMstep em) {
		double[][] marginal;
		marginal = em.convergeExpectation();
		return marginal;
	}
	// --- The Main Method -------------------------------------------------------
	/**
	 * This is the main function where all parameters for the likelihood ratio test can be specified
	 * Please only modify code in this method if you do not plan to change the implementation
	 * Change the inputs according the comments and make sure the working directory is correctly setup
	 * @throws FileNotFoundException 
	 */	
	public static void main ( String[] args ) throws FileNotFoundException {
		
		String currentDir = System.getProperty("user.dir"); //set working directory
		String dir = currentDir + "/data/test/"; //data folder for both input and output
		String input = dir + "test1000g5DC"+".gml"; //input graph, please prepare your data in the GML format
		FileOutputStream output = null;
		String out = dir.concat("test1000Output13DCDC"); //name your output file
		output = new FileOutputStream(out + ".txt");//set output file extensions

		PrintStream print = new PrintStream(output);
		//Tagging the output columns
		MainBayesian test = new MainBayesian(input); //create a MainFunction object for likelihood ratio tests
		print.println("graph #:\t" + "#nodes\t" + "#edges\t\t" + "Group size\t" + "block affinity\t\t\t\t\t\t\t\t\t\t\t\t" + "likelihood"); 
		double[][] typeP = new double[test.graph.getNumType()][test.graph.getNumType()]; //dummy parameters
		double[] gNode = new double[test.graph.getNumType()]; //dummy parameters
		Long start = System.currentTimeMillis(); //timer starts
		for (int k=1; k<13; k++) {
			print.println("groups #:\t" + k); 
			//Set the boolean to true for DC-SBM
			MCMCBayesian bpFast = new MCMCBayesian(test.graph, true, typeP, gNode, 1, k); 
			bpFast.maxLogLike = -Double.MAX_VALUE;; //reset the maximum likelihood
			for (int n=0; n<10; n++) { //set number of test samples
				print.print("graph " + n + ":\t"); //indexing the random graphs for tests
				double[][] marginal = test.posterior(test.graph, bpFast); //call the likeRatio method
				print.print("steps " + marginal[0][0] + ":\t"); 
				print.print("bound " + marginal[1][0] + ":\t"); 
				print.println("ratio:" + bpFast.likelihood); //print the returned log-likelihood difference
			}
		}
		print.println("time:" + (System.currentTimeMillis() - start)); //timer stops and print the time used
		print.close();
	}
}