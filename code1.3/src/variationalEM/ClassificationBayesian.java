package variationalEM;

import graphTools.Graph;

import java.util.Random;

/**
 * This class implements the classification scheme and the measure of likelihood of classifications.
 * Only works for undirected simple graphs with no self-loops
 * 
 * @author Xiaoran Yan
 */

public class ClassificationBayesian {
	
	// --- Instance Variables ----------------------------------------------------
	private int k; //number of groups
	private int[] groups; //group mapping: index-vertex id, value-group id	
	private int[][] aGroup; //adjacent matrix of group connections
	private int[] nGroup; //array for counting vertices in each group
	private int[] dGroup; //array for counting totoal degrees in each group
	static Graph graph; //member graph for edge query
	static double[] logFactTable; //lookup table for factorial calculation
	
	private double[][] groupMatrix; //for p_ij tracking
	public double[][][] avgGroupMatrix;
	
	// --- Constructors ---------------------------------------------------------- 
	/**
	 * This constructor creates a Classification for a given Graph g and the number of groups n. 
	 * When the real flag is true, ground truth labels are used.
	 * @param g Graph
	 * @param n int
	 */
	public ClassificationBayesian(Graph g, int n, boolean real) {
		k = n;
		groups = new int[g.getNumNodes()];
		aGroup = new int[k][k];
		nGroup = new int[k];
		dGroup = new int[k];
		graph = g;
		groupMatrix = new double[k][k];
		
		//Initialize the factorial lookup table
		logFactTable = new double[g.getNumNodes()*g.getNumNodes()+2]; //pending: memory bound
		logFactTable[0] = 0;
		for (int i=1; i<logFactTable.length; i++)
			logFactTable[i] = logFactTable[i-1] + java.lang.Math.log(i);
		int[] topList = new int[groups.length+1];
		if (real) {
			for (int i=0; i<topList.length-1; i++)
				topList[i] = i;
			topList[topList.length-1] = -1;
		}
		else
			topList[0] = -1;
		reInitialize(topList);
	}
	
	/**
	 * This constructor copies a Classification into a new object.
	 * @param c Classification
	 */
	public ClassificationBayesian(ClassificationBayesian c) {
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
	 * This method calculates log factorial. 
	 * @param t int
	 */	
	private double logfact(int t) {
		if (t<=0)
			return 0;
		double f = 0;
		for (int i=1; i<=t; i++)
			f += java.lang.Math.log(i);
		return f;	
	}
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
			groups[list[i]] = graph.vList[list[i]].type % k;
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
			dGroup[m] += graph.vList[i].outDegree;
			//System.out.println("1");
			for (j=0; j<graph.vList[i].targets.size(); j++){
				//System.out.println("4");
				//System.out.println("source-vtxno:  "+i);
				//System.out.println("edges size:  "+graph.vList[i].edges.size());
				//System.out.println("target-vtxno:  "+graph.vList[i].edges.get(j));
				//int target_group=groups[graph.vList[i].edges.get(j)];
				aGroup[m][groups[graph.vList[i].targets.get(j)]]++;
				//if (m==groups[graph.vList[i].targets.get(j)])
					//aGroup[m][m]++;
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
		//for (i=0; i<k; i++)
			//aGroup[i][i] = aGroup[i][i]/2;
		
		//System.out.println("3");
	}
	
	/**
	 * This method returns the likelihood of the current Classification.
	 * Calculation differs according to the properties of the graph (directed, selfloop)
	 * @param null
	 */	
	public double likelihood(double beta, boolean gSize, boolean degreeC) {
		int a, b, i, j;
		double temp = 0;
		groupMatrix = new double[k][k];
		
		//the edge terms for vanilla SBM and DC-SBM
		if (graph.isDirected() && graph.hasSelfloop())
			for (i=0; i<k; i++)
				for (j=0; j<k; j++) {
					a = aGroup[i][j];
					b = nGroup[i] * nGroup[j] - a;
					temp = temp + logFactTable[a] + logFactTable[b] - logFactTable[a+b+1];
					if (a == 0)
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
						temp = temp + logFactTable[a] + logFactTable[b] - logFactTable[a+b+1];
						//temp += logfact(a) + logfact(b) - logfact(a+b-1);
						if (a == 0)
							groupMatrix[i][j] = 0;
						else
							groupMatrix[i][j] = (double) (a/a+b);
					}
					else { //special case for self group
						a = aGroup[i][i] / 2; //avoid duplicate for undirected graphs
						b = nGroup[i] * (nGroup[i]+1) / 2 - a;
						temp = temp + logFactTable[a] + logFactTable[b] - logFactTable[a+b+1];
						//temp += logfact(a) + logfact(b) - logfact(a+b-1);
						if (a == 0)
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
						temp = temp + logFactTable[a] + logFactTable[b] - logFactTable[a+b+1];
						if (a == 0)
							groupMatrix[i][j] = 0;
						else
							groupMatrix[i][j] = (double) a/(a+b);
					}
					else { //special case for self group
						a = aGroup[i][i];
						b = nGroup[i] * (nGroup[i]-1) - a;
						temp = temp + (logFactTable[a] + logFactTable[b] - logFactTable[a+b+1])*3;//factor of 10
						if (a == 0)
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
						temp = temp + logFactTable[a] + logFactTable[b] - logFactTable[a+b+1];
						if (a == 0)
							groupMatrix[i][j] = 0;
						else
							groupMatrix[i][j] = (double) a/(a+b);
					}
					else { //special case for self group
						a = aGroup[i][i] / 2; //avoid duplicate for undirected graphs
						b = nGroup[i] * (nGroup[i]-1) / 2 - a;
						temp = temp + logFactTable[a] + logFactTable[b] - logFactTable[a+b+1];
						if (a == 0)
							groupMatrix[i][j] = 0;
						else
							groupMatrix[i][j] = (double) a/(a+b);
					}
				}
		
		if (degreeC) { //the degree terms for DC-SBM, only works for undirected graphs with self loops. 
			if (!graph.isDirected()) {
				for (i=0; i<k; i++){
					a = nGroup[i]-1;
					if (a < 0)
						a++;
					temp = temp + logFactTable[a] - logFactTable[dGroup[i]+a];
					for (j=0; j<graph.getNumNodes(); j++) {
						if (groups[j] == i)
							temp = temp + logFactTable[graph.vList[j].outDegree];
					}
					temp = temp + (dGroup[i]+nGroup[i]) * Math.log(nGroup[i]);
				}
				temp = temp - Math.log((0.5+graph.getNumNodes()/graph.getNumEdgs()/24)*(graph.getNumNodes()-10)) 
						- 100*Math.log(graph.getNumNodes()); //likelihood normalization with k=10
			}
		}
		if (gSize) { //the group size terms
			temp = temp + logFactTable[k-1] - logFactTable[graph.getNumNodes()+k-1];
			for (i=0; i<k; i++)
				temp += logFactTable[nGroup[i]];
		}
		return beta*temp;
	}
	
	/**
	 * This method returns the distribution of straight likelihood of the heat bath MCMC process,
	 * given a node n and current classification on graph g.
	 * @param n int
	 */	
	public double[] distribution(int n, double beta, boolean gSize, boolean DC) {
			int h, i, j;
			double[] dist = new double[k];
			avgGroupMatrix = new double[k][k][k];
			int orig=groups[n];
			for (i=0; i<k; i++) { //i is the candidate group
				//Classification temp = new Classification(this);				
				mutate(n, i); //new classification if change n into group i
				dist[i] = likelihood(beta, gSize, DC); //get the straight likelihood
				for (j=0; j<k; j++)
					for (h=0; h<k; h++)
						avgGroupMatrix[i][j][h] = avgGroupMatrix[i][j][h] + this.groupMatrix[j][h];
				this.mutate(n, orig);
			}

			double sum = 0; //shrink the numbers proportionally to prevent over(under)flow
			for (i=0; i<k; i++)
				sum = sum + dist[i];
			for (i=0; i<k; i++) {
				dist[i] = dist[i] - sum/(k); //exponents subtract the average 
				if (dist[i]>700) //exponent cut off in case of overflow
					dist[i] = 700;
				dist[i] = java.lang.Math.pow(java.lang.Math.E, dist[i]); //raise it back to likelihood
			}
			sum = 0; //normalize the distribution
			for (i=0; i<k; i++)
				sum = sum +dist[i];
			for (i=0; i<k; i++)
				dist[i] = dist[i] / sum;
			for (i=0; i<k; i++) //get the expected group matrix
				for (j=0; j<k; j++)
					for (h=0; h<k; h++)
						avgGroupMatrix[0][j][h] = avgGroupMatrix[i][j][h] 
						                          + dist[i] * avgGroupMatrix[i][j][h];
			return dist;	
	}
	
	/**
	 * This method changes the vertex n to group c, updating the group adjacency statistics accordingly.
	 * Only works for undirected simple graphs with no self-loops
	 * @param g Graph
	 * @param n int
	 * @param newg int
	 */	
	public void mutate(int n, int newg) {
		int j;
		int oldg = groups[n];
		if (oldg != newg) {//if the group has changed

			for (j = 0; j < graph.vList[n].sources.size(); j++) {
				int type = groups[graph.vList[n].sources.get(j)];
				aGroup[type][oldg]--;
				//if (type != oldg)
					aGroup[oldg][type]--;
				aGroup[type][newg]++;
				//if (type != newg)
					aGroup[newg][type]++;;
			}
			
			groups[n] = newg; // set the new group
			nGroup[oldg]--; // adjust the vertex count
			dGroup[oldg] -= graph.vList[n].outDegree;// update the total degree count (only works for undirected graphs)
			nGroup[newg]++;			
			dGroup[newg] += graph.vList[n].outDegree;
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