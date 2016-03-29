package variationalEM;

import graphTools.Graph;

import java.util.ArrayList;
import java.util.Random;


/**
 * This class implements the linear BP with mean-field non-edge messages
 * A child class extends the abstract EMstep class
 * Built for undirected multi-graphs and Poisson/DC block models
 * 
 * @author Xiaoran Yan ( everyxt@gmail.com )
 * @version BP_1.3
 * @time Nov, 2013
 */

public class BPfastMU extends EMstep{
	//--- Inner class for messages -----------------------------------------------
	public static class Vmsg { //message for one node
		public ArrayList<double[]> msgTargets; //message following and reversing directed edges
		public double[] others; //message to none-connected nodes, also serves as the marginal vectors for each node
		
		private Vmsg() { //inner constructor
			msgTargets = new ArrayList<double[]>();
		}
	}
	// --- Instance Variables ----------------------------------------------------
	Vmsg[] message; //message passing on each vertex along edges / non-edges
	// --- Constructors ----------------------------------------------------------
	public BPfastMU(){} //the empty constructor
	/**
	 * This constructor creates a linear BP with random initial messages and manually set block parameters
	 * @param g Graph
	 * @param degreeC boolean
	 * @param gSize boolean
	 * @param typeP double[][]
	 * @param gNode double[]
	 */
	public BPfastMU(Graph g, boolean degreeC, boolean gSize,  double[][] typeP, double[] gNode) {
		super(g, degreeC, gSize, typeP, gNode); //calling the parent constructor		
		epsilon = 0.000001 * graph.getNumEdgs() * graph.getNumType(); //overwrites default threshold

		Random r = new Random(); //random initialization of messages
		message = new Vmsg[graph.getNumNodes()];
		for (int i=0; i<message.length; i++) {
			message[i] = new Vmsg();
			//initializing none-edge messages
			message[i].others = new double[graph.getNumType()];
			double sum=0;
			for (int k=0; k<graph.getNumType(); k++) {
				double rand = r.nextDouble();						
				message[i].others[k] = (5+rand); //balanced initialization
				sum += message[i].others[k];
			}
			for (int k=0; k<graph.getNumType(); k++)
				message[i].others[k] /= sum; //normalization
			
			//initializing directed messages
			int msgs = graph.vList[i].targets.size() + graph.vList[i].sources.size();
			for (int j=0; j<msgs; j++){
				double[] msg = new double[graph.getNumType()];
				sum=0;
				for (int k=0; k<msg.length; k++) {
					double rand = r.nextDouble();						
					msg[k] = (5+rand); //balanced initialization
					sum += msg[k];
				}
				for (int k=0; k<msg.length; k++)
					msg[k] /= sum; //normalization
				message[i].msgTargets.add(msg);
			}
		}
	}
	
	/**
	  * This constructor creates a linear BP from an exact copy, including block parameters as well as messages
	 * @param g Graph
	 * @param copy BPfastMU
	 */
	public BPfastMU(Graph g, BPfastMU copy) {
		super(g, copy.degreeCorrect, copy.gSizeCorrect, copy.typeP, copy.gNode); //calling the parent constructor

		likelihood = copy.likelihood;
		likelihoodHard = copy.likelihoodHard;
		message = new Vmsg[graph.getNumNodes()];
		for (int i=0; i<message.length; i++) {
			message[i] = new Vmsg();
			//copying none-edge messages / marginal vectors
			message[i].others = copy.message[i].others;
			//copying directed messages
			int msgsN = graph.vList[i].targets.size() + graph.vList[i].sources.size();
			for (int j=0; j<msgsN; j++){
				double[] msg = new double[graph.getNumType()];
				for (int k=0; k<msg.length; k++) 			
					msg[k] = copy.message[i].msgTargets.get(j)[k];
				message[i].msgTargets.add(msg);
			}
		}
		
	}
	
	// --- Instance Methods ------------------------------------------------------
	/**
	 * This method does a sweep of message updates (asynchronous) across the network, mean-fielding the non-edge messages
	 * @return measure of change in terms of messages in L1 norm
	 * @param null
	 */
	public double stepBPfast() {
		//base message for speed up (case 2: no edges)
		double[] baseMsg = new double[graph.getNumType()];
		
		if (degreeCorrect) { //for the DC model
			for (int h=0; h<graph.getNumNodes(); h++) {
				for (int k1=0; k1<graph.getNumType(); k1++) {
					double temp = 0;
					for (int k2=0; k2<graph.getNumType(); k2++) {
						temp += message[h].others[k2] * typeP[k2][k1];
					}
					baseMsg[k1] += temp * (graph.vList[h].outDegree); //note that since we have bi-directed edges for undirected graphs, outDegree = degree
				}
			}
		}
		else { //for the vanilla model
			for (int h=0; h<graph.getNumNodes(); h++) {
				for (int k1=0; k1<graph.getNumType(); k1++) {
					double temp = 0;
					for (int k2=0; k2<graph.getNumType(); k2++) {
						temp += message[h].others[k2] * java.lang.Math.exp(-typeP[k2][k1]);
					}
					baseMsg[k1] += java.lang.Math.log(temp);
				}
			}
		}

		int[] Ulist = permute(); //pick a random update order
		double delta = 0; //measure of change
		// Update in order according to the permutation
		for (int i=0; i<message.length; i++)
			for (int j=0; j<message[Ulist[i]].msgTargets.size()+1; j++) { //plus 1 for the non-edge messages
				int totalD = 1;
				if (degreeCorrect) //for the DC model
					totalD = graph.vList[Ulist[i]].outDegree; //note that since we have bi-directed edges for undirected graphs, outDegree = degree
				double[] newmsg = new double[graph.getNumType()];
				for (int k1=0; k1<newmsg.length; k1++) {//updating message from i to j
					if (degreeCorrect) //for the DC model
						newmsg[k1] = -baseMsg[k1] * totalD;
					else
						newmsg[k1] = baseMsg[k1];
				}
				
				for (int k1=0; k1<graph.getNumType(); k1++) { //get rid of duplicate self term
					double temp = 0;
					for (int k2=0; k2<graph.getNumType(); k2++) {
						temp += message[Ulist[i]].others[k2] 
						    * java.lang.Math.exp(-typeP[k2][k1] * totalD * totalD);
					}
					if (temp != 0)
						newmsg[k1] = newmsg[k1] - java.lang.Math.log(temp);
				}
				
				int target = -1; //index of the message target
				if (j < message[Ulist[i]].msgTargets.size()) { //get rid of duplicate target term
					if (j <graph.vList[Ulist[i]].targets.size())
						target = graph.vList[Ulist[i]].targets.get(j);
					else //the message target is on a reversed edge
						target = graph.vList[Ulist[i]].sources.get(j-graph.vList[Ulist[i]].targets.size());
					for (int k1=0; k1<graph.getNumType(); k1++) {
						double temp = 0;
						for (int k2=0; k2<graph.getNumType(); k2++) {
							temp += message[target].others[k2] 
							    * java.lang.Math.exp(-typeP[k2][k1] * (graph.vList[target].outDegree) * totalD);
						}
						if (temp != 0)
							newmsg[k1] = newmsg[k1] - java.lang.Math.log(temp);
					}
				}
				//Ready to taking account for the directed messages on observed edges
				for (int l=0; l<graph.vList[Ulist[i]].sources.size(); l++) { //neighboring message following the edge
					double temp1 = 0;
					double temp2 = 0;
					int source = graph.vList[Ulist[i]].sources.get(l); //source node from the edge
					int index = graph.vList[source].targets.indexOf(Ulist[i]); //index in the msgTargets[source]
					int totalD2 = 1; //default vanilla model
					if (degreeCorrect) //for the DC model
						totalD2 = graph.vList[source].outDegree;
					if (j < message[Ulist[i]].msgTargets.size()) { // (i,j) in E
						if (source!=Ulist[i] && source!=target) { //avoid self and target messages
							if (graph.vList[Ulist[i]].targets.contains(source)) { //case 1: double edges between Ulist[i] and source	
								for (int k1=0; k1<graph.getNumType(); k1++) {
									temp1 = 0;
									temp2 = 0;
									for (int k2=0; k2<graph.getNumType(); k2++) {
										int edgeC = graph.vList[source].targetCount.get(index);
										temp1 += message[source].msgTargets.get(index)[k2] 
										    * poisson(typeP[k2][k1] * totalD2 * totalD, edgeC)
											* java.lang.Math.exp(-typeP[k2][k1] * totalD * totalD2);
										temp2 += message[source].others[k2] 
										    * java.lang.Math.exp(-typeP[k2][k1] * totalD * totalD2);
									}
									if (temp2 != 0)
										newmsg[k1] = newmsg[k1] - java.lang.Math.log(temp2) + java.lang.Math.log(temp1);
								}
							} //case 2: no edges, do nothing. It is already handled by the base message by construction.
						}
					}
					else { // (i,j) not in E
						if (source != Ulist[i]) { //avoid self message
							if (graph.vList[Ulist[i]].targets.contains(source)) { //case 1: double edges between Ulist[i] and source	
								for (int k1=0; k1<graph.getNumType(); k1++) {
									temp1 = 0;
									temp2 = 0;
									for (int k2=0; k2<graph.getNumType(); k2++) {
										int edgeC = graph.vList[source].targetCount.get(index);
										temp1 += message[source].msgTargets.get(index)[k2] 
										    * poisson(typeP[k2][k1] * totalD2 * totalD, edgeC)
											* java.lang.Math.exp(-typeP[k2][k1] * totalD * totalD2);
										temp2 += message[source].others[k2] 
										    * java.lang.Math.exp(-typeP[k2][k1] * totalD * totalD2);
										}
									if (temp2 != 0)
										newmsg[k1] = newmsg[k1] - java.lang.Math.log(temp2) + java.lang.Math.log(temp1);
								}
							} //case 2: no edges, do nothing. It is already handled by the base message by construction.
						}
					}
				}
				
				double sum = 0;
				for (int k=0; k<graph.getNumType(); k++) {
					if (gSizeCorrect) //group size correction
						newmsg[k] = java.lang.Math.exp(newmsg[k] + java.lang.Math.log(gNode[k]));
					else
						newmsg[k] = java.lang.Math.exp(newmsg[k]);
					sum += newmsg[k]; //for message normalization
				}
				for (int k=0; k<graph.getNumType(); k++){
					if (sum == 0) //boundary cases
						newmsg[k] = 1.0/graph.getNumType();
					else
						newmsg[k] = newmsg[k] / sum;
					if (j < message[Ulist[i]].msgTargets.size()) {
						delta += java.lang.Math.abs(message[Ulist[i]].msgTargets.get(j)[k]-newmsg[k]);
						message[Ulist[i]].msgTargets.get(j)[k] = 0.5*message[Ulist[i]].msgTargets.get(j)[k] + 0.5*newmsg[k]; //damping propagation
					}
					else {
						delta += java.lang.Math.abs(message[Ulist[i]].others[k]-newmsg[k]);
						message[Ulist[i]].others[k] = 0.5*message[Ulist[i]].others[k] + 0.5*newmsg[k]; //damping propagation
					}
				}
				
			}
		return delta;
	}
	
	/**
	  * This method implements the E-step inner-loops for the linear BP
	  * A polymorphic extension of the abstract method convergeExpectation in parent class
	  * Built for undirected multi-graphs and Poisson/DC block models
	  * @return block marginal vectors for all nodes
	  * @param null
	 */
	public double[][] convergeExpectation() {
		//Initializations
		steps = 0;
		double delta = Double.MAX_VALUE;
		//The inner loop with a bound of 10 sweeps
		while(delta > epsilon && steps <10) {
			delta = stepBPfast(); //Do a E-step sweep across all nodes
			steps++;
		}
		double[][] marginal = new double[graph.getNumNodes()][graph.getNumType()];
		for (int i=0; i<marginal.length; i++) // non-edge message = marginal vector for each node
			for (int k=0; k<graph.getNumType(); k++) 
				marginal[i][k] = message[i].others[k];
		return marginal;
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
		likelihood = 0; //initialization
		double likeTemp = 0;
		double[][] temp = new double[graph.getNumType()][graph.getNumType()];
		double[][] sum = new double[graph.getNumType()][graph.getNumType()];	
		
		for (int i=0; i<graph.getNumNodes(); i++) {
			int iid = findMax(message[i].others); //for hard block assignment
			double[][] q = new double[graph.getNumType()][graph.getNumType()]; //pair-wise marginals
			int totalDi = 1; //default vanilla model
			if (degreeCorrect) //for DC model
				totalDi = graph.vList[i].outDegree;
			double totalDj = 1.0; //default vanilla model
			double[] gNodeNE = new double[graph.getNumType()]; //for tracking average point-wise non-edge marginals
			
			for (int j=0; j<graph.vList[i].targets.size()+1; j++) {	//plus 1 for the non-edge messages
				int edgeC = 0; //default no edge
				if (j<graph.vList[i].targets.size()) { //case 1: i->j in E
					int jIndex = graph.vList[i].targets.get(j); //pick the target
					int jid = findMax(message[jIndex].others); //for hard block assignment	 			
					if (degreeCorrect) //for DC model
						totalDj = graph.vList[jIndex].outDegree;
					int source = graph.vList[jIndex].sources.indexOf(i)+graph.vList[jIndex].targets.size(); //index in the message list
					edgeC = graph.vList[i].targetCount.get(j); //get edge count
					
					for (int k1=0; k1<graph.getNumType(); k1++) {
						for (int k2=0; k2<graph.getNumType(); k2++) { //product to get pair-wise marginal				
							q[k1][k2] =  message[i].msgTargets.get(j)[k1] * message[jIndex].msgTargets.get(source)[k2]
								    * poisson(typeP[k1][k2] * totalDi * totalDj, edgeC)
									* Math.exp(- typeP[k1][k2] * totalDi * totalDj);
							if (jIndex>i && k1 == iid && k2 == jid) //for hard block assignment	 		
								likelihoodHard += Math.log (poisson(typeP[k1][k2] * totalDi * totalDj, edgeC))
									+(- typeP[k1][k2] * totalDi * totalDj);
						}				
						gNodeNE[k1] +=  message[j].others[k1]; //sum to get point-wise marginals
					}
				}
				else { //case 2: no edge (Mean field approximation applied on all non-edge pairs)
					if (degreeCorrect) //use average degree for DC model over all non-edge pairs
						totalDj = graph.getNumEdgs()*2.0 / graph.getNumNodes();
					
					for (int k1=0; k1<graph.getNumType(); k1++) {
						double jMsgAvgk1 = marginals[k1] * (graph.getNumNodes()-1) - gNodeNE[k1]; //mean-field approximation on non-edge
						for (int k2=0; k2<graph.getNumType(); k2++) {							
							q[k2][k1] = message[i].others[k2] * jMsgAvgk1 / (graph.getNumNodes()-j-1)
							        * Math.exp(- typeP[k2][k1] * totalDi * totalDj);
							if (k2 == iid) //for hard block assignment	
								likelihoodHard +=  jMsgAvgk1 *0.5* (- typeP[k2][k1] * totalDi * totalDj);
						}
					}
				}
				
				double qSum = 0; //q normalization			
				for (int k1=0; k1<graph.getNumType(); k1++)
					for (int k2=0; k2<graph.getNumType(); k2++) 
						qSum += q[k1][k2];					
				for (int k1=0; k1<graph.getNumType(); k1++)
					for (int k2=0; k2<graph.getNumType(); k2++) 
						q[k1][k2] = q[k1][k2] / qSum;
				
				for (int k1=0; k1<graph.getNumType(); k1++) //accumulating sums for estimating typeP
					for (int k2=0; k2<graph.getNumType(); k2++) {
						double pLike = 0; //pairwise part
						int NedgeC = 1; //counting non-edge pairs
						
						if (j<graph.vList[i].targets.size()) { //case 1: i->j in E
							int jIndex = graph.vList[i].targets.get(j);
							if (jIndex>i) {
								pLike += q[k1][k2]* (Math.log(poisson(typeP[k1][k2] * totalDi * totalDj, edgeC))
													   -typeP[k1][k2] * totalDi * totalDj); //pairwise energy part
								pLike -= q[k1][k2]* Math.log(q[k1][k2]); //pairwise entropy part
							}
						}
						else { //case 1: i->j not in E
							NedgeC = (graph.getNumNodes()-1-j);
							pLike += NedgeC *0.5* q[k1][k2]* (-typeP[k1][k2] * totalDi * totalDj); //pairwise energy part
							pLike -= q[k1][k2]* Math.log(q[k1][k2]); //pairwise entropy part
						}
						temp[k1][k2] += NedgeC * q[k1][k2]* edgeC;
						sum[k1][k2] += NedgeC * q[k1][k2] * totalDi * totalDj;
							
						likeTemp += pLike;
						//}
					}
			}
			double pEntr = 0; //pointwise part
			for (int k=0; k<graph.getNumType(); k++) {				
				pEntr += message[i].others[k] * Math.log(message[i].others[k]); //pointwise entropy part
				if (gSizeCorrect) {
					likeTemp += message[i].others[k] * Math.log(gNode[k]); //pointwise energy part (group size correction)
					if (k == iid)
						likelihoodHard += Math.log(gNode[k]); //group size correction for hard block assignment
				}
			}
			likeTemp += pEntr*(graph.vList[i].targets.size()-1); //pointwise entropy scaling (+1 for the dummy none-edge node)
		}
		double pEntr2 = 0; //pointwise entropy for the dummy none-edge node
		for (int k=0; k<graph.getNumType(); k++)
			pEntr2 += gNode[k] * Math.log(gNode[k]);
		likeTemp += pEntr2 * graph.getNumNodes();
		likelihood = likeTemp;
		
		double delta = 0; //change measure
		if (!fix) //do not update block parameters if M-step is fixed
			for (int k1=0; k1<graph.getNumType(); k1++) {
				gNode[k1] = marginals[k1];
				for (int k2=0; k2<graph.getNumType(); k2++) {
					if (sum[k1][k2] == 0)
						temp[k1][k2] = 0;
					else // normalization
						temp[k1][k2] = temp[k1][k2] / sum[k1][k2];
					delta += Math.abs(temp[k1][k2]-typeP[k1][k2]);
					//Handling boundary cases
					if (temp[k1][k2] < 1e-128)
						typeP[k1][k2] = 1e-128;
					else if (temp[k1][k2] > (1 - 1e-16))
						typeP[k1][k2] = 1-1e-16;
					else
						typeP[k1][k2] = temp[k1][k2];
				}
			}
		return delta;
	}

}