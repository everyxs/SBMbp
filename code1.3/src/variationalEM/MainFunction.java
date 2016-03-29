package variationalEM;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.Random;
import graphTools.Graph;

/**
 * This class includes the parameters and initializations for using the variational EM framework 
 * for likelihood ratio test between the vanilla block model and degree corrected (DC) block model.
 * This example showcases how different models can be fitted to the same data sets under this framework, and
 * how the polymorphism mechanism can be used to employ different EM implementations.
 * Please only edit the main function at the bottom for testing purposes. 
 * For other uses of the variational EM framework, users can write their own main functions according the API.
 * 
 * @author Xiaoran Yan ( everyxt@gmail.com )
 * @version BP_1.3
 * @time Nov, 2013
 */
public class MainFunction {
	// --- Instance Variables ----------------------------------------------------
	public Graph graph; //input graph
	public EMstep innerMAP; //the maximum likelihood EM setup (abstract class to be instantiated)
	public EMiterate outerMAP; //the maximum likelihood EM setup (outer-loop)
	public double[][] typeP; //edge block parameters for tuning purposes
	double[] gNode; //node block parameters for tuning purposes
	public double likelihood; //the maximum likelihood
	public double likelihoodVanilla; //the maximum likelihood for vanilla model
	
	// --- Constructors ---------------------------------------------------------- 
	public MainFunction() {} //the empty constructor
	
	/**
	 * This constructor initializes internal parameters given the input
	 * @param dir String
	 */
	public MainFunction(String dir) {
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
		typeP = new double[graph.getNumType()][graph.getNumType()];
		gNode = new double[graph.getNumType()];
		
		int[] typeNodeNum = graph.getNumNodesType();
		int[][] typeEdgeNum = graph.getNumEdgsType();
		double[][] typeDegreeNum = graph.getNumDegreeType();
		// Ground truth values for the block parameters
		for (int i=0; i<typeP.length; i++)
			for (int j=0; j<typeP[i].length; j++) {
				if (graph.degreeCorrect)
					typeP[i][j] = (double)typeDegreeNum[i][j] / typeNodeNum[i] / typeNodeNum[j];
				else
					typeP[i][j] = (double)typeEdgeNum[i][j] / typeNodeNum[i] / typeNodeNum[j];
			}
		for (int i=0; i<gNode.length; i++){
			gNode[i] = (double) typeNodeNum[i] / graph.getNumNodes();
		}
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
		int[] typeNodeNum = newGraph.getNumNodesType();
		int[][] typeEdgeNum = newGraph.getNumEdgsType();
		double[][] typeDegreeNum = newGraph.getNumDegreeType();
		for (int i=0; i<typeP.length; i++)
			for (int j=0; j<typeP[i].length; j++) {
				typeP[i][j] = (double)typeEdgeNum[i][j] / typeNodeNum[i] / typeNodeNum[j];
			}
		for (int i=0; i<gNode.length; i++){
			gNode[i] = (double) typeNodeNum[i] / newGraph.getNumNodes();
		}
		restart(newGraph, em); //grid search parameters values for maximum likelihood
		//manual(newGraph, em); //start around the ground truth parameters
		//tune(newGraph, em); //search based on grid search results 
		fineTune(newGraph, em); //fine-grained search using minute perturbations
		likelihoodVanilla = likelihood; //remember the maximum likelihood
		for (int i=0; i<newGraph.getNumType(); i++) { //print the final parameters
			print.print("["+innerMAP.gNode[i]+"]"+"\t");
			for (int j=0; j<newGraph.getNumType(); j++) {
				print.print(innerMAP.typeP[i][j]+"\t");
				print.print(outerMAP.typeP[i][j]+"\t");					
			}
		}
		print.print(likelihood+"\t");
		
		// Fitting the degree corrected (DC) block model
		em.degreeCorrect = true;
		print.print("degree corrected" + "\t");
		typeNodeNum = newGraph.getNumNodesType();
		typeEdgeNum = newGraph.getNumEdgsType();
		typeDegreeNum = newGraph.getNumDegreeType();
		for (int i=0; i<typeP.length; i++) //initializing parameters around the ground truth parameters (unused for grid search)
			for (int j=0; j<typeP[i].length; j++) {
				typeP[i][j] = (double)typeDegreeNum[i][j] / typeNodeNum[i] / typeNodeNum[j];
			}
		for (int i=0; i<gNode.length; i++){
			gNode[i] = (double) typeNodeNum[i] / newGraph.getNumNodes();
		}
		restart(newGraph, em); //grid search parameters values for maximum likelihood
		//manual(newGraph, em); //start around the ground truth parameters
		//tune(newGraph, em); //search based on grid search results 
		fineTune(newGraph, em); //fine-grained search using minute perturbations
		for (int i=0; i<graph.getNumType(); i++) { //print the final parameters
			print.print("["+innerMAP.gNode[i]+"]"+"\t");
			for (int j=0; j<graph.getNumType(); j++) {
				print.print(innerMAP.typeP[i][j]+"\t");
				print.print(outerMAP.typeP[i][j]+"\t");		
			}
		}
		print.print(likelihood+"\t");
		return likelihood - likelihoodVanilla; //return the log-likelihood difference
	}
	/**
	 * This method passes manually specified parameters values to the EM framework and 
	 * calls the convergeEMfix method which does NOT change parameters throughout the EM iterations
	 * Use this method for testing purpose
	 * @param g Graph
	 * @param em EMstep (abstract class to be instantiated)
	 */
	public void manualFix(Graph g, EMstep em) {
		EMiterate outer = new EMiterate(em, 0.1, 0.1); //instantiate an EMiterate according to em

		//See the initialization in the constructors
		double like = outer.convergeEMfix(gNode, typeP); //passes the parameters
		likelihood = like; //update the most likely trackers
		innerMAP = em;
		outerMAP = new EMiterate(innerMAP, outer);
	}
	/**
	 * This method passes manually specified parameters values to the EM framework and 
	 * calls the convergeEM method which does change parameters throughout the EM iterations
	 * Use this method when there is enough domain knowledge
	 * @param g Graph
	 * @param em EMstep (abstract class to be instantiated)
	 */
	public void manual(Graph g, EMstep em) {
		EMiterate outer = new EMiterate(em, 0.1, 0.1); //instantiate an EMiterate according to em

		//See the initialization in the constructors
		double like = outer.convergeEM(gNode, typeP); //passes the parameters
		likelihood = like; //update the most likely trackers
		innerMAP = em;
		outerMAP = new EMiterate(innerMAP, outer);
	}
	/**
	 * This method does a grid search for parameters values based on the EM framework and 
	 * calls the convergeEMgrid method which requires no specified parameter inputs
	 * Use this method when there is no domain knowledge about the block structures
	 * @param g Graph
	 * @param em EMstep (abstract class to be instantiated)
	 */
	public void restart(Graph g, EMstep em) {
		
		double[] assort = new double[4]; //grid search coordinates
		assort[0] = 0.0005;
		assort[1] = 0.005;
		assort[2] = 0.05;
		assort[3] = 0.3;
		//assort[4] = 0.6;
		likelihood = -Double.MAX_VALUE;
		
		// The grid search
		for (int i=0; i<assort.length; i++) { //try different initial assortativeness
			for (int j=0; j<assort.length; j++) { //try different initial dis-assortativeness
				for (int k=0; k<3; k++) { //multiple runs
					System.out.println("restart "+i+j);
					EMiterate outer = new EMiterate(em, assort[i], assort[j]); //instantiate an EMiterate according to em
					double like = outer.convergeEMgrid();
					if (like >= likelihood) { //higher likelihood found and update the trackers
						innerMAP = em;
						outerMAP = new EMiterate(innerMAP, outer);
						likelihood = like;
					}
				}
			}
		}
	}
	/**
	 * This method changes both node and edge block parameters and explores the space beyond the most likely result
	 * calls the convergeEM method which does change parameters throughout the EM iterations
	 * Use this method after restart or manual(Fix) for parameter space exploration
	 * @param g Graph
	 * @param em EMstep (abstract class to be instantiated)
	 */
	public void tune(Graph g, EMstep em) {
		
		for (int i=0; i<10; i++) { //number of attempts
			double[][] P = new double[g.getNumType()][g.getNumType()];
			double[] M = new double[g.getNumType()];
			P = exploreP(innerMAP.typeP, true); //targeted new edge parameters
			M = exploreM(innerMAP.gNode); //targeted new note parameters
			EMiterate outer = new EMiterate(em, 0.1, 0.1); //instantiate an EMiterate according to em
			double like = outer.convergeEM(M, P); //pass the targeted parameters
			if (like > likelihood) { //higher likelihood found and update the trackers
				innerMAP = em;
				outerMAP = new EMiterate(innerMAP, outer);
				likelihood = like;
			}
		}
	}
	/**
	 * This method changes only edge block parameters and explores the space around the most likely result
	 * calls the convergeEM method which does change parameters throughout the EM iterations
	 * Use this method last for fine-grained parameter space exploration around the current maximum
	 * @param g Graph
	 * @param em EMstep (abstract class to be instantiated)
	 */
	public void fineTune(Graph g, EMstep em) {
		
		for (int i=0; i<5; i++) { //number of attempts
			double[][] P = new double[g.getNumType()][g.getNumType()];
			P = exploreP(innerMAP.typeP, false); //targeted new edge parameters (fine-grained)
			EMiterate outer = new EMiterate(em, 0.1, 0.1); //instantiate an EMiterate according to em
			double like = outer.convergeEM(em.gNode, P); //pass the targeted edge parameters with existing node parameters
			if (like > likelihood) { //higher likelihood found and update the trackers
				innerMAP = em;
				outerMAP = new EMiterate(innerMAP, outer);
				likelihood = like;
			}
		}
	}
	/**
	 * This is the method that explores the node parameter space randomly
	 * Used in the method tune
	 * @return targeted new note parameters double[]
	 * @param M double[]
	 */
	public double[] exploreM(double[] M) {
		double[] newM = new double[M.length]; //targeted new note parameters
		for (int i=0; i<M.length; i++)
			newM[i] = M[i]; //use the original note parameters as default
		Random r = new Random();
		double sum=0;
		for (int i=0; i<M.length; i++) { //overwrite with random values 
			newM[i] = r.nextDouble();
			sum += newM[i]; 
		}
		for (int i=0; i<M.length; i++)
			newM[i] = newM[i]/sum; //re-normalization
		return newM;
	}
	/**
	 * This is the method that explores the edge parameter space randomly or in close neighborhood
	 * Used in the method tune and fineTune
	 * @return targeted new edge parameters double[][]
	 * @param M double[]
	 */	
	public double[][] exploreP(double[][] P, boolean rand) {
		double[] assort = new double[6];
		assort[0] = 0.00001; //grid search coordinates
		assort[1] = 0.00025;
		assort[2] = 0.005;
		assort[3] = 0.05;
		assort[4] = 0.2;
		assort[5] = 0.5;
		//assort[6] = 0.6;
		//assort[7] = 1;
		
		double[][] newP = new double[P.length][P.length]; //targeted new edge parameters
		for (int i=0; i<P.length; i++)
			for(int j=0; j<P.length; j++)
				newP[i][j] = P[i][j]; //use the original edge parameters as default
		Random r = new Random();
		if (rand) //if we are doing the random exploration
			for (int i=0; i<innerMAP.graph.getNumType(); i++) { //do k restarts
				int a = r.nextInt(P.length);
				int b = r.nextInt(P.length);
				newP[a][b] = assort[r.nextInt(assort.length)]; //overwrite with random values 
			}
		else //if we are doing fine-grained search
			for (int i=0; i<innerMAP.graph.getNumType(); i++) { //do k restarts
				int a = r.nextInt(P.length);
				int b = r.nextInt(P.length);
				if (r.nextDouble()<0.5)
					newP[a][b] = P[a][b]*0.9; //make a minute perturbation to the current MLE
				else {
					newP[a][b] = P[a][b]*1.1; //make a minute perturbation to the current MLE
					if (newP[a][b] > 1) //boundary cases
						newP[a][b] = 1;
				}
			}
		return newP;
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
		String input = dir + "test1000g2"+".gml"; //input graph, please prepare your data in the GML format
		FileOutputStream output = null;
		String out = dir.concat("test1000Output13"); //name your output file
		output = new FileOutputStream(out + ".txt");//set output file extensions

		PrintStream print = new PrintStream(output);
		//Tagging the output columns
		print.println("graph #:\t" + "#nodes\t" + "#edges\t\t" + "Group size\t" + "block affinity\t\t\t\t\t\t\t\t\t\t\t\t" + "likelihood"); 
		
		MainFunction test = new MainFunction(input); //create a MainFunction object for likelihood ratio tests
		Long start = System.currentTimeMillis(); //timer starts
		for (int n=0; n<10; n++) { //set number of test samples
			print.print("graph " + n + ":\t"); //indexing the random graphs for tests
			EMstep bpFast = new BPfastMU(test.graph, false, true, test.typeP, test.gNode); //choose a specific EM algorithm by instantiating a child class
			double ratio = test.likeRatio(print, bpFast); //call the likeRatio method
			print.print("ratio:" + ratio); //print the returned log-likelihood difference
			print.println();
		}
		print.println("time:" + (System.currentTimeMillis() - start)); //timer stops and print the time used
	}
}