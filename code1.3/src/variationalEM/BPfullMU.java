package variationalEM;

import graphTools.Graph;
import java.util.Random;


/**
 * This class implements the full BP with all pair messages (slower but more accurate result for smaller graphs)
 * A child class extends the abstract EMstep class
 * Built for undirected multi-graphs and Poisson/DC block models
 * 
 * @author Xiaoran Yan ( everyxt@gmail.com )
 * @version BP_1.3
 * @time Nov, 2013
 */

public class BPfullMU extends EMstep{
	
	// --- Instance Variables ----------------------------------------------------
	double[][][] message; //message passing on each edge / none-edge
	double[][] marginal; //the mixed membership vector, correspond to non-edge messages in BPfast
	// --- Constructors ----------------------------------------------------------
	public BPfullMU() {} //the empty constructor
	/**
	 * This constructor creates a linear BP with random initial messages and random block parameters
	 * @param g Graph
	 * @param degreeC boolean
	 * @param gSize boolean
	 */
	public BPfullMU(Graph g, boolean degreeC, boolean gSize) {
		super(g, degreeC, gSize);
		Random r = new Random();
		
		message = new double[graph.getNumNodes()][graph.getNumNodes()][graph.getNumType()];
		marginal = new double[graph.getNumNodes()][graph.getNumType()];
		for (int i=0; i<message.length; i++)
			for (int j=0; j<message[0].length; j++){
				double sum=0;
				for (int k=0; k<message[0][0].length; k++) {
					double rand = r.nextDouble();						
					message[i][j][k] = 2+rand;
					sum += message[i][j][k];
				}
				for (int k=0; k<message[0][0].length; k++)
					message[i][j][k] = message[i][j][k]/sum;
			}
		//message[0][33][0] = 10;
	}
	/**
	 * This constructor creates a linear BP with random initial messages and manually set block parameters
	 * @param g Graph
	 * @param degreeC boolean
	 * @param gSize boolean
	 * @param typeP double[][]
	 * @param gNode double[]
	 */
	public BPfullMU(Graph g, boolean degree, boolean gSize,  double[][] typeP, double[] n) {
		super(g, degree, gSize, typeP, n);
		Random r = new Random();
		
		epsilon = 0.0001 * graph.getNumEdgs() * graph.getNumType();

		likelihoodHard=0;
		like3=0;
		message = new double[graph.getNumNodes()][graph.getNumNodes()][graph.getNumType()];
		marginal = new double[graph.getNumNodes()][graph.getNumType()];
		for (int i=0; i<message.length; i++)
			for (int j=0; j<message[0].length; j++){
				double sum=0;
				for (int k=0; k<message[0][0].length; k++) {
					double rand = r.nextDouble();						
					message[i][j][k] = 2+rand;
					sum += message[i][j][k];
				}
				for (int k=0; k<message[0][0].length; k++)
					message[i][j][k] = message[i][j][k]/sum;
			}
	}
	
	/**
	  * This constructor creates a full BP from an exact copy, including block parameters as well as messages
	 * @param g Graph
	 * @param copy BPfastMU
	 */
	public BPfullMU(Graph g, BPfullMU copy) {
		super(g, copy.degreeCorrect, copy.gSizeCorrect, copy.typeP, copy.gNode);
		
		likelihoodHard = copy.likelihoodHard;
		like3 = copy.like3;
		for (int i=0; i<gNode.length; i++)
			gNode[i] = copy.gNode[i];
		message = new double[graph.getNumNodes()][graph.getNumNodes()][graph.getNumType()];
		marginal = new double[graph.getNumNodes()][graph.getNumType()];

		for (int i=0; i<message.length; i++)
			for (int j=0; j<message[0].length; j++)
				for (int k=0; k<message[0][0].length; k++)
					message[i][j][k] = copy.message[i][j][k];
		for (int j=0; j<message[0].length; j++)
			for (int k=0; k<message[0][0].length; k++)
				marginal[j][k] = copy.marginal[j][k];
	}

	// --- Instance Methods ------------------------------------------------------
	/**
	 * This method does a sweep of all message updates (asynchronous) across the network
	 * @return measure of change in terms of messages in L1 norm
	 * @param null
	 */
	public double stepBP() {
		int[] Ulist = permute();
		double delta = 0;
		
		// update in order according to the permutation
		for (int i=0; i<message.length; i++) {
			int totalD = 1;
			if (degreeCorrect)
				totalD = graph.vList[Ulist[i]].outDegree;
			for (int j=0; j<message[0].length; j++) if (j!= Ulist[i]) {
				int totalD2 = 1;
				if (degreeCorrect)
					totalD2 = graph.vList[j].outDegree;
				double[] newmsg = new double[message[0][0].length];
				for (int k=0; k<newmsg.length; k++) //updating message from i to j
					newmsg[k] = 1;
				
				for (int l=0; l<message[0].length; l++)	if (l!=Ulist[i] && l!=j) {
					double temp = 0;
					// loop through all the other incoming neighbors
					if (graph.vList[Ulist[i]].sources.contains(l))	{//case 1: double edges between Ulist[i] and l
						int index = graph.vList[l].targets.indexOf(Ulist[i]); //index in the msgTargets[source]
						int edgeC = graph.vList[l].targetCount.get(index);
						for (int k1=0; k1<message[0][0].length; k1++) {
							temp = 0;
							for (int k2=0; k2<message[0][0].length; k2++) {
								temp += message[l][Ulist[i]][k2] 
									    * poisson(typeP[k2][k1] * totalD2 * totalD, edgeC)
										* java.lang.Math.exp(-typeP[k2][k1] * totalD * totalD2);
							
							}
							newmsg[k1] *= temp;
						}
					}
					else {//case 4: no edge between Ulist[i] and the l
						for (int k1=0; k1<message[0][0].length; k1++) {
							temp = 0;
							for (int k2=0; k2<message[0][0].length; k2++) {
									temp += message[l][Ulist[i]][k2] //degree correction
											* java.lang.Math.exp(-typeP[k2][k1] * totalD * totalD2);
								}
							newmsg[k1] *= temp;
						}
					}
				}
				double sum = 0;
				for (int k=0; k<message[0][0].length; k++) {
					if (gSizeCorrect) //group size correction
						newmsg[k] = newmsg[k]*gNode[k];
					sum += newmsg[k]; //for message normalization
				}
				for (int k=0; k<message[0][0].length; k++){
					if (sum == 0)
						newmsg[k] = 1.0/message[0][0].length;
					else
						newmsg[k] = newmsg[k] / sum;
					delta += java.lang.Math.abs(message[Ulist[i]][j][k]-newmsg[k]);
					message[Ulist[i]][j][k] = 0.5*message[Ulist[i]][j][k] + 0.5*newmsg[k]; //damping propagation
				}
			}
		}
		return delta;
	}
	
	/**
	 * This method does a sweep of all message updates (synchronous flood) across the network
	 * @return measure of change in terms of messages in L1 norm
	 * @param null
	 */
	public double stepBPflood() {
		// Newly updated message after a single flood step
		double[][][] newmsg = new double[message.length][message[0].length][message[0][0].length];
		for (int i=0; i<message.length; i++)
			for (int j=0; j<message[0].length; j++)
				for (int k=0; k<message[0][0].length; k++)
					newmsg[i][j][k] = 1;
		// Message flood
		for (int i=0; i<message.length; i++) {
			int totalD = 1;
			if (degreeCorrect)
				totalD = graph.vList[i].outDegree;
			for (int j=0; j<message[0].length; j++) if (j!=i){
					int totalD2 = 1;
					if (degreeCorrect)
						totalD2 = graph.vList[j].outDegree;
				for (int l=0; l<message[0].length; l++)	if (l!=i && l!=j){
					int edgeC = 0;
					if (graph.vList[l].targets.contains(i)) {
						int target = graph.vList[l].targets.indexOf(i);
						edgeC = graph.vList[l].targetCount.get(target);	
					}
					double temp = 0;
					// loop through all the other incoming neighbors
					for (int k1=0; k1<message[0][0].length; k1++) {
						temp = 0;
						for (int k2=0; k2<message[0][0].length; k2++) {
							temp += message[l][i][k2] 
								    * poisson(typeP[k2][k1] * totalD2 * totalD, edgeC)
									* java.lang.Math.exp(-typeP[k2][k1] * totalD * totalD2);
						
						}
						newmsg[i][j][k1] *= temp;
					}
				}
			}
		}
		
		// Update all the messages and measures the change
		double delta = 0;
		for (int i=0; i<message.length; i++) 
			for (int j=0; j<message[0].length; j++) if (j!=i) {
				double sum = 0;
				for (int k=0; k<message[0][0].length; k++) {
					if (gSizeCorrect)
						newmsg[i][j][k] = newmsg[i][j][k]*gNode[k]; //group size correction
					sum += newmsg[i][j][k]; // for message normalization
				}
				for (int k=0; k<message[0][0].length; k++){
					if (sum == 0)
						newmsg[i][j][k] = 1.0/message[0][0].length;
					else
						newmsg[i][j][k] = newmsg[i][j][k] / sum;
					delta += java.lang.Math.abs(message[i][j][k]-newmsg[i][j][k]); 
					message[i][j][k] = 0.5*message[i][j][k] + 0.5*newmsg[i][j][k]; //damping propagation
				}
			}
		return delta;
	}
	
	/**
	  * This method implements the M-step after the E-step converges
	  * A polymorphic extension of the abstract method mStep in the parent class
	  * Built for undirected multi-graphs and Poisson/DC block models
	  * @return maximum likelihood after this EM iteration
	  * @param marginals double[]
	  * @param fix boolean
	 */
	public double mStep(double[] n, boolean fix) {
		likelihood = 0;
		likelihoodHard = 0;
		like3 = 0;
		double[][] temp = new double[graph.getNumType()][graph.getNumType()];
		double[][] sum = new double[graph.getNumType()][graph.getNumType()];
		for (int i=0; i<graph.getNumNodes(); i++) {
			int iid = findMax(marginal[i]);
			for (int j=i+1; j<graph.getNumNodes(); j++) {
				int jid = findMax(marginal[j]);
				int edgeC = 0;
				if (graph.vList[i].targets.contains(j)) {
					int target = graph.vList[i].targets.indexOf(j);
					edgeC = graph.vList[i].targetCount.get(target);	
				}
				int totalDi = 1;
				if (degreeCorrect)
					totalDi = graph.vList[i].outDegree;
				int totalDj = 1;
				if (degreeCorrect)
					totalDj = graph.vList[j].outDegree;
				double[][] q = new double[graph.getNumType()][graph.getNumType()];
				for (int k1=0; k1<graph.getNumType(); k1++)
					for (int k2=0; k2<graph.getNumType(); k2++) {
						q[k1][k2] = message[i][j][k1] * message[j][i][k2] * poisson(typeP[k1][k2] * totalDi * totalDj, edgeC)
									* java.lang.Math.exp(-typeP[k1][k2] * totalDi * totalDj);
						if (k1 == iid && k2 == jid)
							like3 += java.lang.Math.log (poisson(typeP[k1][k2] * totalDi * totalDj, edgeC)
									* java.lang.Math.exp(-typeP[k1][k2] * totalDi * totalDj));
					}
				double qSum = 0; //q normalization
			
				for (int k1=0; k1<graph.getNumType(); k1++)
					for (int k2=0; k2<graph.getNumType(); k2++) 
						qSum += q[k1][k2];
				
				for (int k1=0; k1<graph.getNumType(); k1++)
					for (int k2=0; k2<graph.getNumType(); k2++) {
						q[k1][k2] = q[k1][k2] / qSum;
					}
				//likelihood += java.lang.Math.log(qSum);	
				
				for (int k1=0; k1<graph.getNumType(); k1++) //accumulating sums for estimating typeP
					for (int k2=0; k2<graph.getNumType(); k2++) {
						//if (k1!=k2||i!=j) {
						double pLike = 0; //pairwise part
						pLike -= q[k1][k2]*java.lang.Math.log(q[k1][k2]); //pairwise entropy part
							
						pLike += q[k1][k2]*java.lang.Math.log(poisson(typeP[k1][k2] * totalDi * totalDj, edgeC)
															* java.lang.Math.exp(-typeP[k1][k2] * totalDi * totalDj)); //pairwise energy part
						temp[k1][k2] += q[k1][k2]* edgeC;
						sum[k1][k2] += q[k1][k2] * totalDi * totalDj;
							
						likelihoodHard += pLike;
						//}
					}
			}
			double pEntr = 0; //pointwise part
			for (int k=0; k<graph.getNumType(); k++) {
				pEntr += marginal[i][k] * java.lang.Math.log(marginal[i][k]); //pointwise entropy part
				if (gSizeCorrect) {
					likelihoodHard += marginal[i][k]*(java.lang.Math.log(gNode[k])); //pointwise energy part
					if (k == iid)
						like3 += java.lang.Math.log(gNode[k]); //group size correction
				}
			}
			likelihoodHard += pEntr*(graph.getNumNodes()-2);
		}

		double delta = 0; //change measure
		if (!fix)
			for (int k1=0; k1<graph.getNumType(); k1++) {
				gNode[k1] = n[k1];
				for (int k2=0; k2<graph.getNumType(); k2++) {
					if (sum[k1][k2] == 0)
						temp[k1][k2] = 0;
					else // normalization
						temp[k1][k2] = temp[k1][k2] / sum[k1][k2];
					delta += java.lang.Math.abs(temp[k1][k2]-typeP[k1][k2]);
					if (temp[k1][k2] < 1e-128)
						typeP[k1][k2] = 1e-128;
					else if (temp[k1][k2] > (1 - 1e-16))
						typeP[k1][k2] = 1-1e-16;
					else 
						typeP[k1][k2] = temp[k1][k2];
				}
			}
		likelihood = likelihoodHard;
		return delta;
	}
	
	/**
	  * This method implements the E-step inner-loops for the full BP
	  * A polymorphic extension of the abstract method convergeExpectation in parent class
	  * Built for undirected multi-graphs and Poisson/DC block models
	  * @return block marginal vectors for all nodes
	  * @param null
	 */
	public double[][] convergeExpectation() {
		steps = 0;
		double delta = Double.MAX_VALUE;
		while(delta > epsilon && steps <10) {
			delta = stepBPflood();
			steps++;
		}
		double temp = 0;
		for (int i=0; i<message.length; i++) {
			int j = (i+1)%message.length; // pick an arbitrary pair
			int totalDi = 1;
			if (degreeCorrect)
				totalDi = graph.vList[i].outDegree;
			int totalDj = 1;
			if (degreeCorrect)
				totalDj = graph.vList[j].outDegree;
			int edgeC = 0;
			if (graph.vList[i].targets.contains(j)) {
				int target = graph.vList[i].targets.indexOf(j);
				edgeC = graph.vList[i].targetCount.get(target);	
			}
			for (int k1=0; k1<message[0][0].length; k1++) {
				temp = 0;
				for (int k2=0; k2<message[0][0].length; k2++)
					temp += message[j][i][k2]* poisson(typeP[k2][k1] * totalDi * totalDj, edgeC)
								* java.lang.Math.exp(-typeP[k2][k1] * totalDi * totalDj);
				marginal[i][k1] = temp*message[i][j][k1];
			}
		}
		//Normalization
		for (int i=0; i<message.length; i++) {
			double sum = 0;
			for (int k=0; k<message[0][0].length; k++)
				sum += marginal[i][k];
			for (int k=0; k<message[0][0].length; k++) {
				if (sum == 0)
					marginal[i][k] = 1.0/marginal[0].length;
				else
					marginal[i][k] = marginal[i][k]/sum;
			}
		}
		return marginal;
	}
	
}