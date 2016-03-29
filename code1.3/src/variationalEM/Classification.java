package variationalEM;

import graphTools.Graph;

import java.util.Random;

/**
 * This class implements the classification scheme and the measure of likelihood of classifications.
 * Subclass used under the MAP MCMC Estep
 * Built for undirected multi-graphs and Poisson/DC block models
 * 
 * @author Xiaoran Yan ( everyxt@gmail.com )
 * @version BP_1.3
 * @time Nov, 2013
 */

public class Classification {
	
	// --- Instance Variables ----------------------------------------------------
	private int k; //number of groups
	public int[] groups; //group mapping: index-vertex id, value-group id	
	public int[][] aGroup; //adjacent matrix of group connections
	public int[] nGroup; //array for counting vertices in each group
	private int[] dGroup; //array for counting totoal degrees in each group
	static Graph graph; //member graph for edge query
	
	public double[][] groupMatrix; //for p_ij update
	public double[][][] avgGroupMatrix;
	
	// --- Constructors ---------------------------------------------------------- 
	/**
	 * This constructor creates a Classification for a given Graph g and the number of groups n.
	 * @param g Graph
	 * @param n int
	 */
	public Classification(Graph g, int n) {
		k = n;
		groups = new int[g.getNumNodes()];
		aGroup = new int[k][k];
		nGroup = new int[k];
		dGroup = new int[k];
		graph = g;
		groupMatrix = new double[k][k];

		int[] topList = new int[1];
		topList[0] = -1;
		reInitialize(topList);
	}
	
	/**
	 * This constructor copies a Classification into a new object.
	 * @param c Classification
	 */
	public Classification(Classification c) {
		int i, j;
		k = c.k;
		groups = new int[c.groups.length];
		aGroup = new int[k][k];
		nGroup = new int[k];
		dGroup = new int[k];
		groupMatrix = new double[k][k];
		
		for (i=0; i<c.groups.length; i++)
			groups[i] = c.groups[i];
		for (i=0; i<k; i++) {
			nGroup[i] = c.nGroup[i];
			dGroup[i] = c.dGroup[i];
			for (j=0; j<k; j++)
				aGroup[i][j] = c.aGroup[i][j];
		}
		int[] topList = new int[1];
		topList[0] = -1;
		reInitialize(topList);
	}
	
	// --- Instance Methods ------------------------------------------------------	
	/**
	 * This method randomly initialize a Classification for a given Graph.
	 * @param list int[]
	 */
	public void reInitialize(int[] list) {
		int i, j, m;
		//Reinitialize random group assignment 
		//System.out.println("1");
		Random r = new Random();
		for (i=0; i<graph.getNumNodes(); i++)
			groups[i] = r.nextInt(k);
		i = 0;
		//System.out.println("2");
		while (list[i] != -1) { //overwrite with the true classification
			groups[list[i]] = graph.vList[list[i]].type;
			i++;
		}
		//System.out.println("3");
		//Reinitialize group statistics
		for (i=0; i<k; i++) {
			nGroup[i] = 0;
			dGroup[i] = 0;
			for (j=0; j<k; j++)
				aGroup[i][j] = 0;
		}
		//System.out.println("k: "+k);
		for (i=0; i<graph.getNumNodes(); i++) {
			m = groups[i];
			nGroup[m]++;
			dGroup[m] += graph.vList[i].outDegree + graph.vList[i].inDegree;
			//System.out.println("1");
			for (j=0; j<graph.vList[i].targets.size(); j++){
				//System.out.println("4");
				//System.out.println("source-vtxno:  "+i);
				//System.out.println("edges size:  "+graph.vList[i].edges.size());
				//System.out.println("target-vtxno:  "+graph.vList[i].edges.get(j));
				//int target_group=groups[graph.vList[i].edges.get(j)];
				aGroup[m][groups[graph.vList[i].targets.get(j)]] += graph.vList[i].targetCount.get(j);
				//if(m>5)
					//System.out.println("m: "+m);
				//if(groups[graph.vList[i].edges.get(j)]>5)
					//System.out.println("mm: "+groups[graph.vList[i].edges.get(j)]);
				//if(graph.vList[i].edges.get(j)>487)
					//System.out.println("target: "+graph.vList[i].edges.get(j));
				//System.out.println("2");
			}
			//System.out.println("3");
		}
		//System.out.println("3");
	}
	
	/**
	 * This method returns the log-likelihood of the current Classification.
	 * Calculation differs according to the properties of the graph (directed, selfloop)
	 * @param beta double
	 * @param gSize boolean
	 * @param degreeC boolean
	 * @param typeP double[][]
	 * @param gNode double[]
	 */	
	public double likelihood(double beta, boolean gSize, boolean degreeC, double[][] typeP, double[] gNode) {
		int a, b, i, j;
		double temp = 0;
		groupMatrix = new double[k][k];
		
		for (i=0; i<typeP.length; i++)
			for (j=0; j<typeP[i].length; j++) {
				if (typeP[i][j] < 1e-128)
					typeP[i][j] = 1e-128;
				else if(typeP[i][j] > (1 - 1e-16))
					typeP[i][j] = 1 - 1e-16;
			}
		temp = 0;		
		for (i=0; i<gNode.length; i++){
			if (gNode[i] < 1e-128)
				gNode[i] = 1e-128;
			else if(gNode[i] > (1 - 1e-16))
				gNode[i] = 1 - 1e-16;
		}
		
		temp = 0;
		if (!degreeC) {
			if (graph.isDirected() && graph.hasSelfloop())
				for (i=0; i<k; i++)
					for (j=0; j<k; j++) {
						a = aGroup[i][j];
						b = nGroup[i] * nGroup[j] - a;
						temp = temp + a*java.lang.Math.log(typeP[i][j]) + b*java.lang.Math.log(1-typeP[i][j]);
						if (a+b == 0)
							groupMatrix[i][j] = 0;
						else
							groupMatrix[i][j] = (double) a/(a+b);
					}
			if (!graph.isDirected() && graph.hasSelfloop())
				for (i=0; i<k; i++)
					for (j=i; j<k; j++) { //avoid duplicate for undirected graphs
						if (i != j) {
							a = aGroup[i][j];
							b = nGroup[i] * nGroup[j] - a;
							temp = temp + a*java.lang.Math.log(typeP[i][j]) + b*java.lang.Math.log(1-typeP[i][j]);
							if (a+b == 0)
								groupMatrix[i][j] = 0;
							else
								groupMatrix[i][j] = (double) (a/a+b);
						}
						else { //special case for self group
							a = aGroup[i][i] / 2; //avoid duplicate for undirected graphs
							b = nGroup[i] * (nGroup[i]+1) / 2 - a;
							temp = temp + a*java.lang.Math.log(typeP[i][j]) + b*java.lang.Math.log(1-typeP[i][j]);
							if (a+b == 0)
								groupMatrix[i][j] = 0;
							else
								groupMatrix[i][j] = (double) a/(a+b);
						}
					}
			if(graph.isDirected() && !graph.hasSelfloop())
				for (i=0; i<k; i++)
					for (j=0; j<k; j++) {
						if (i != j) {
							a = aGroup[i][j];
							b = nGroup[i] * nGroup[j] - a;
							temp = temp + a*java.lang.Math.log(typeP[i][j]) + b*java.lang.Math.log(1-typeP[i][j]);
							if (a+b == 0)
								groupMatrix[i][j] = 0;
							else
								groupMatrix[i][j] = (double) a/(a+b);
						}
						else { //special case for self group
							a = aGroup[i][i];
							b = nGroup[i] * (nGroup[i]-1) - a;
							temp = temp + a*java.lang.Math.log(typeP[i][j]) + b*java.lang.Math.log(1-typeP[i][j]);
							if (a+b == 0)
								groupMatrix[i][j] = 0;
							else
								groupMatrix[i][j] = (double) a/(a+b);
						}
					}
			if(!graph.isDirected() && !graph.hasSelfloop())
				for (i=0; i<k; i++)
					for (j=i; j<k; j++) { //avoid duplicate for undirected graphs
						if (i != j) {
							a = aGroup[i][j];
							b = nGroup[i] * nGroup[j] - a;
							temp = temp + a*java.lang.Math.log(typeP[i][j]) + b*java.lang.Math.log(1-typeP[i][j]);
							if (a+b == 0)
								groupMatrix[i][j] = 0;
							else
								groupMatrix[i][j] = (double) a/(a+b);
						}
						else { //special case for self group
							a = aGroup[i][i] / 2; //avoid duplicate for undirected graphs
							b = nGroup[i] * (nGroup[i]-1) / 2 - a;
							temp = temp + a*java.lang.Math.log(typeP[i][j]) + b*java.lang.Math.log(1-typeP[i][j]);
							if (a+b == 0)
								groupMatrix[i][j] = 0;
							else
								groupMatrix[i][j] = (double) a/(a+b);
						}
					}
		}
		/**
		 *  the multi-graph DC likelihood only works for directed (now) graphs with self loops. 
		 */
		else {
			if (graph.isDirected()) {
				for (i=0; i<k; i++)
					for (j=0; j<k; j++) {
						a = aGroup[i][j];
						b = dGroup[i] * dGroup[j];
						temp += a * java.lang.Math.log((double) a/b);
						if (a+b == 0)
							groupMatrix[i][j] = 0;
						else
							groupMatrix[i][j] = (double) a/b;
					}
			}
		}
		if (gSize)
			for (i=0; i<k; i++)
				temp += nGroup[i]*java.lang.Math.log(gNode[i]);
		return beta*temp;
	}
	
	/**
	 * This method returns the distribution of straight likelihood of the heat bath MCMC process,
	 * given a node n and current classification on graph g.
	 * @param n int
	 * @param beta double
	 * @param gSize boolean
	 * @param typeP double[][]
	 * @param gNode double[]
	 */	
	public double[] distribution(int n, double beta, boolean gSize, boolean DC, double[][] typeP, double[] gNode) {
			int h, i, j;
			double[] dist = new double[k];
			avgGroupMatrix = new double[k][k][k];
			int orig=groups[n];
			for (i=0; i<k; i++) { //i is the candidate group
				//Classification temp = new Classification(this);				
				mutate(n, i); //new classification if change n into group i
				dist[i] = likelihood(beta, gSize, DC, typeP, gNode); //get the straight log-likelihood
				for (j=0; j<k; j++)
					for (h=0; h<k; h++)
						avgGroupMatrix[i][j][h] = this.groupMatrix[j][h];
				this.mutate(n, orig);
			}

			double sum = 0; //shrink the numbers proportionally to prevent over(under)flow
			for (i=0; i<k; i++)
				sum = sum + dist[i];
			for (i=0; i<k; i++) {
				dist[i] = dist[i] - sum/(k); //exponents subtract the average 
				if (dist[i]>700) //exponent cut off in case of overflow
					dist[i] = 700;
				dist[i] = java.lang.Math.exp(dist[i]); //raise it back to likelihood
			}
			sum = 0; //normalize the distribution
			for (i=0; i<k; i++)
				sum = sum +dist[i];
			for (i=0; i<k; i++)
				dist[i] = dist[i] / sum;
			
			for (j=0; j<k; j++) //get the expected group matrix
				for (h=0; h<k; h++)
					avgGroupMatrix[0][j][h] = dist[0] * avgGroupMatrix[0][j][h];
			for (i=1; i<k; i++)
				for (j=0; j<k; j++)
					for (h=0; h<k; h++)
						avgGroupMatrix[0][j][h] += dist[i] * avgGroupMatrix[i][j][h];
			return dist;	
	}
	
	/**
	 * This method changes the vertex n to group c, updating the group adjacency statistics accordingly.
	 * @param n int
	 * @param newg int
	 */	
	public void mutate(int n, int newg) {
		int j;
		int oldg = groups[n];
		if (oldg != newg) {//if the group has changed
			//for edges to other nodes
			for (j = 0; j < graph.vList[n].sources.size(); j++) {
				int type = groups[graph.vList[n].sources.get(j)];
				aGroup[type][oldg] -= graph.vList[n].sourceCount.get(j);
				aGroup[type][newg] += graph.vList[n].sourceCount.get(j);
			}
			for (j = 0; j < graph.vList[n].targets.size(); j++) {
				int type = groups[graph.vList[n].targets.get(j)];
				aGroup[oldg][type] -= graph.vList[n].targetCount.get(j);
				aGroup[newg][type] += graph.vList[n].targetCount.get(j);
			}// for self loops (needs update)
			if (graph.hasSelfloop()) {
				if (graph.vtxSelfloopFlag[n] == 1) {
					aGroup[newg][newg]++;
					aGroup[oldg][newg]--;
					aGroup[oldg][oldg]++;
					aGroup[newg][oldg]--;
				}
			}
	
			groups[n] = newg; // set the new group
			nGroup[oldg]--; // adjust the vertex count
			nGroup[newg]++;			
		}
	}
	
	/**
	 * This method returns the number of vertices.
	 * @param null
	 */	
	public int numGroup() { return k; }
	
	/**
	 * This method returns the type of vertex.
	 * @param vtxno int
	 */			
	public int getVtxType(int vtxno){ return groups[vtxno];	}	
	
	public void printmodel() {
		System.out.println("vertices:");
		for(int i=0;i<graph.getNumNodes();i++){
			System.out.println("vertex "+i+": "+groups[i]);
		}
		//System.out.println("groups:");
		System.out.println("group edges:");
		for(int i=0;i<k;i++){
			for(int j=0;j<k;j++){
				System.out.println(i+"  and  "+j+":  "+aGroup[i][j]);
			}
		}
	}
}