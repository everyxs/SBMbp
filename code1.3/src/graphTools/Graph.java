package graphTools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;
import java.util.StringTokenizer;


/**
 * This class implements the adjacency list representation of graphs. 
 * 
 * The parser in this version assumes input files have correct GML format,
 * assuming no duplicated VERTEX declarations (duplicated edges allowed);
 * 
 * @author Xiaoran Yan & Yaojia Zhu ( everyan@cs.unm.edu, yaojia.zhu@gmail.com )
 * @version Likelihood_mutual_Information_MCMC_1.0
 * @time May 18, 2009
 */

public class Graph {
	
	//--- Inner class for vertices -----------------------------------------------
	public class Vertex {
		//private String name;
		public String value;
		public int type;		
		public String id;
		public int index;		
		public ArrayList<Integer> targets;
		public ArrayList<Integer> targetCount;
		public int outDegree;
		public ArrayList<Integer> sources;
		public ArrayList<Integer> sourceCount;
		public int inDegree;
		
		private Vertex() {
			id = "";
			value = "";	
			targets = new ArrayList<Integer>();
			targetCount = new ArrayList<Integer>();
			outDegree = 0;
			sources = new ArrayList<Integer>();
			sourceCount = new ArrayList<Integer>();
			inDegree = 0;
		}
		public boolean addTarget(int target) {
			if(targets.contains(target)) {
				if (multiEdge) {
					int index = targets.indexOf(target);
					int countOld = targetCount.get(index).intValue();
					targetCount.set(index, countOld+1);
					outDegree++;
					return true;
				}
				else return false;
			}
			else {
				targets.add(target);
				targetCount.add(1);
				outDegree++;
				return true;
			}
		}
		public boolean addSource(int source) {
			if(sources.contains(source)) {
				if (multiEdge) {
					int index = sources.indexOf(source);
					int countOld = sourceCount.get(index).intValue();
					sourceCount.set(index, countOld+1);
					inDegree++;
					return true;
				}
				else return false;
			}
			else {
				sources.add(source);
				sourceCount.add(1);
				inDegree++;
				return true;
			}
		}
	}
	
	// --- Instance Variables ----------------------------------------------------
	public Vertex[] vList; //list of vertices
	public int[] vtxSelfloopCount; //indicate which vertices have self-loop 
	public ArrayList<String> listIndex2Id = new ArrayList<String>();
	public ArrayList<String> listType2Value = new ArrayList<String>();
	public ArrayList<Integer> listIndex2Type = new ArrayList<Integer>();
	public ArrayList<Double> listIndex2Value2 = new ArrayList<Double>();
	public int numEgs, numVtx, numType, numTrueType;
	public double[][] typeP;
	public double[] gNode;
	private boolean directed;
	private boolean selfLoop;
	public boolean multiEdge;
	public boolean degreeCorrect;
	
	// --- Constructors ---------------------------------------------------------- 
	public Graph (){}
	/**
	 * This constructor creates an graph from an input file of GML format.
	 * @param input FileReader
	 */
	public Graph (FileReader input, boolean DC) { //additional parameters required to manipulate group number
		//numType = types;
		listIndex2Id = new ArrayList<String>();
		listType2Value = new ArrayList<String>();
		listIndex2Type = new ArrayList<Integer>();		
		try {
			//Use buffering, reading one line at a time
			BufferedReader buffer =  new BufferedReader(input);
		    try {
		    	String word="";
		    	String line="";
		    	String sid="";
		    	String svalue="";
		    	directed = false;
		    	selfLoop = false; //set the selfLoop parameter
		    	multiEdge = false;
		    	degreeCorrect = DC;
		    	StringTokenizer tokens;
		    	numVtx = 0; //number of vertices
		    	numEgs = 0; //number of edges
		    	numTrueType =0; //number of types
		    	int source=0, target=0;
		    	boolean end = false;
		    	boolean invalidedge = false;
		        
		    	while ((line = buffer.readLine()) != null) { //read one line
		    		//System.out.println("line: "+line);
		    		//Break the line into words
		    		tokens = new StringTokenizer(line);
		    		if (tokens.hasMoreTokens())
		    			word = tokens.nextToken();
		    		//Check if the graph is directed
		    		if (word.equals("directed"))
		    			if (tokens.nextToken().equals("1"))
		    				directed = true;
		    		//Count the number of vertices
		    		if (word.equals("id")){
		    			sid=line.trim().substring(2).trim();
		    			//System.out.println("sid: "+sid);
		    			if(!(listIndex2Id.contains(sid)))
		    				listIndex2Id.add(sid);
		    			else{
		    				while(!(line = buffer.readLine()).trim().equals("]"))
		    					;
		    			}
		    		}
		    		//Read the value field (true classification)
		    		if (word.equals("value")){
		    			svalue=line.trim().substring(5).trim();	
		    			//System.out.println("svalue: "+svalue);
		    			if(listType2Value.contains(svalue)){
		    				int itype=listType2Value.indexOf(svalue);
		    				listIndex2Type.add(itype);
		    			}
		    			else{
		    				listType2Value.add(svalue);		    				
		    				listIndex2Type.add(numTrueType++);
		    			}		    			
		    		    numVtx++;
		    		}
		    		
		    		if (word.equals("TrophicLevel")){
		    			//System.out.println("svalue: "+svalue);
	    				listIndex2Value2.add(Double.valueOf(tokens.nextToken()));	
		    		}
		    		
		    		//Finish counting, create the array of vertices
		    		if (word.equals("edge") && end == false) {
		    			numType = numTrueType;
		    			if (end == false) {
		    				vtxSelfloopCount=new int[numVtx];
		    				vList = new Vertex[numVtx]; 
		    				System.out.println("numVtx: "+numVtx);
		    				for (int i=0; i<numVtx; i++) {
		    					vtxSelfloopCount[i]=0;
		    					vList[i] = new Vertex(); //initialize the edge lists
		    					vList[i].index = i;
		    					vList[i].id = listIndex2Id.get(i);
		    					vList[i].type=listIndex2Type.get(i);
		    					vList[i].value=listType2Value.get(vList[i].type);
		    				}
		    			}		    		
		    			end = true;
		    		}
		    		//Get the source
		    		if (word.equals("source")) {
		    			if (tokens.hasMoreTokens())
		    				sid = tokens.nextToken();	
		    			source = listIndex2Id.indexOf(sid);
		    				//System.out.println("source index: "+source);
		    			if(source==-1)
		    				invalidedge=true;
		    			else 
		    				invalidedge=false;		    				
		    		}
		    		//Add the edge from source to target
		    		if (word.equals("target")) {
		    			if (tokens.hasMoreTokens())
		    				sid = tokens.nextToken();
		    			target = listIndex2Id.indexOf(sid);
		    			//System.out.println("target index: "+target);
		    			if(invalidedge||(target==-1))
		    				continue;
		    			if(source==target){ //eliminate self loops
		    				vtxSelfloopCount[source]++;
		    				//numEgs++;
		    				selfLoop = true;
		    			}
			    		else {
			    			if (vList[source].addTarget(target)) {
				    			vList[target].addSource(source);
				    			numEgs++; //number of edges		
			    			}
			    		}
		    			if (!directed){ //undirected -> bi-directed
		    				if (vList[target].addTarget(source)) {
			    				vList[source].addSource(target);
			    				//numEgs++;
		    				}
		    			}
		    		}
		    	}
		    }
		    finally {
		    	input.close();
		    	System.out.println("Num of Nodes: "+getNumNodes());
		    	System.out.println("Num of Type: "+getNumType());
		    	System.out.println("Num of Edges: "+getNumEdgs());
		    	//
				/*System.out.println("sources and targets:");
				for(int i=0;i<this.getNumNodes();i++){
					System.out.println("vertex "+i+" sources: ");
					for (int j=0; j<this.vList[i].sources.size(); j++){
						System.out.println(this.vList[i].sources.get(j).toString());
					}
					//System.out.println();
					System.out.println("vertex "+i+" targets: ");
					for (int j=0; j<this.vList[i].targets.size(); j++){
						System.out.println(this.vList[i].targets.get(j).toString());
					}
					//System.out.println();			
				}	*/	    	
		    }
		}
		catch (IOException ex) {
			ex.printStackTrace();
		}
		
		
		typeP = new double[getNumType()][getNumType()];
		gNode = new double[getNumType()];
		int[] typeNodeNum = getNumNodesType();
		int[][] typeEdgeNum = getNumEdgsType();
		double[][] typeDegreeNum = getNumDegreeType();
		for (int i=0; i<typeP.length; i++)
			for (int j=0; j<typeP[i].length; j++) {
				//if (degreeCorrect)
					//typeP[i][j] = (double)typeDegreeNum[i][j] / typeNodeNum[i] / typeNodeNum[j];
				//else
					typeP[i][j] = (double)typeEdgeNum[i][j] / typeNodeNum[i] / typeNodeNum[j];
			}
		
		for (int i=0; i<gNode.length; i++){
			gNode[i] = (double) typeNodeNum[i] / getNumNodes();
		}
	}
	
	/**
	 * This constructor creates an graph with only nodes.
	 * @param input FileReader
	 */
	public Graph (FileReader input) { //additional parameters required to manipulate group numbership
		//numType = types;
		listIndex2Id = new ArrayList<String>();
		listType2Value = new ArrayList<String>();
		listIndex2Type = new ArrayList<Integer>();		
		try {
			//Use buffering, reading one line at a time
			BufferedReader buffer =  new BufferedReader(input);
		    try {
		    	String word="";
		    	String line="";
		    	String sid="";
		    	String svalue="";
		    	directed = false;
		    	selfLoop = false; //set the selfLoop parameter
		    	multiEdge = false;
		    	degreeCorrect = false;
		    	StringTokenizer tokens;
		    	numVtx = 0; //number of vertices
		    	numEgs = 0; //number of edges
		    	numTrueType =0; //number of types
		    	boolean end = false;
		        
		    	while ((line = buffer.readLine()) != null) { //read one line
		    		//System.out.println("line: "+line);
		    		//Break the line into words
		    		tokens = new StringTokenizer(line);
		    		if (tokens.hasMoreTokens())
		    			word = tokens.nextToken();
		    		//Check if the graph is directed
		    		if (word.equals("directed"))
		    			if (tokens.nextToken().equals("1"))
		    				directed = true;
		    		//Count the number of vertices
		    		if (word.equals("id")){
		    			sid=line.trim().substring(2).trim();
		    			//System.out.println("sid: "+sid);
		    			if(!(listIndex2Id.contains(sid)))
		    				listIndex2Id.add(sid);
		    			else{
		    				
		    			}
		    		}
		    		//Read the value field (true classification)
		    		if (word.equals("value")){
		    			svalue=line.trim().substring(5).trim();	
		    			//System.out.println("svalue: "+svalue);
		    			if(listType2Value.contains(svalue)){
		    				int itype=listType2Value.indexOf(svalue);
		    				listIndex2Type.add(itype);
		    			}
		    			else{
		    				listType2Value.add(svalue);		    				
		    				listIndex2Type.add(numTrueType++);
		    			}		    			
		    		    numVtx++;
		    		}
		    		
		    		if (word.equals("TrophicLevel")){
		    			//System.out.println("svalue: "+svalue);
	    				listIndex2Value2.add(Double.valueOf(tokens.nextToken()));	
		    		}
		    		
		    		//Finish counting, create the array of vertices
		    		if (word.equals("edge") && end == false) {
		    			numType = numTrueType;
		    			if (end == false) {
		    				vtxSelfloopCount=new int[numVtx];
		    				vList = new Vertex[numVtx]; 
		    				System.out.println("numVtxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx: "+numVtx);
		    				for (int i=0; i<numVtx; i++) {
		    					vtxSelfloopCount[i]=0;
		    					vList[i] = new Vertex(); //initialize the edge lists
		    					vList[i].index = i;
		    					vList[i].id = listIndex2Id.get(i);
		    					vList[i].type=listIndex2Type.get(i);
		    					vList[i].value=listType2Value.get(vList[i].type);
		    				}
		    			}		    		
		    			end = true;
		    		}
		    	}
		    }
		    finally {
		    	input.close();
		    	System.out.println("Num of Type: "+getNumType());
		    	System.out.println("Num of Edges: "+getNumEdgs());
		    	//
				/*System.out.println("sources and targets:");
				for(int i=0;i<this.getNumNodes();i++){
					System.out.println("vertex "+i+" sources: ");
					for (int j=0; j<this.vList[i].sources.size(); j++){
						System.out.println(this.vList[i].sources.get(j).toString());
					}
					//System.out.println();
					System.out.println("vertex "+i+" targets: ");
					for (int j=0; j<this.vList[i].targets.size(); j++){
						System.out.println(this.vList[i].targets.get(j).toString());
					}
					//System.out.println();			
				}	*/	    	
		    }
		}
		catch (IOException ex) {
			ex.printStackTrace();
		}
		
		typeP = new double[getNumType()][getNumType()];
		gNode = new double[getNumType()];
		int[] typeNodeNum = getNumNodesType();
		int[][] typeEdgeNum = getNumEdgsType();
		double[][] typeDegreeNum = getNumDegreeType();
		for (int i=0; i<typeP.length; i++)
			for (int j=0; j<typeP[i].length; j++) {
				//if (degreeCorrect)
					//typeP[i][j] = (double)typeDegreeNum[i][j] / typeNodeNum[i] / typeNodeNum[j];
				//else
					typeP[i][j] = (double)typeEdgeNum[i][j] / typeNodeNum[i] / typeNodeNum[j];
			}
		
		for (int i=0; i<gNode.length; i++){
			gNode[i] = (double) typeNodeNum[i] / getNumNodes();
		}
	}
	
	/**
	 * This constructor creates an exact graph copy.
	 * @param input FileReader
	 */
	public Graph (Graph input) { //additional parameters required to manipulate group numbership
		listIndex2Id = input.listIndex2Id;
		listType2Value = input.listType2Value;
		listIndex2Type = input.listIndex2Type;
		listIndex2Value2 = input.listIndex2Value2;
		numEgs = input.numEgs;
		numVtx = input.numVtx;
		numType = input.numType;
		numTrueType = input.numTrueType;
		directed = input.directed;
		selfLoop = input.selfLoop;
		multiEdge = input.multiEdge;
		degreeCorrect = input.degreeCorrect;
		typeP = input.typeP;
		gNode = input.gNode;
		
		vtxSelfloopCount=new int[numVtx];
		vList = new Vertex[numVtx];
		for (int i=0; i<numVtx; i++) {
			vtxSelfloopCount[i]=input.vtxSelfloopCount[i];
			vList[i] = new Vertex(); //initialize the edge lists
			vList[i].index = input.vList[i].index;
			vList[i].id = input.vList[i].id;
			vList[i].type = input.vList[i].type;
			vList[i].value = input.vList[i].value;
			for (int j=0; j<input.vList[i].targets.size(); j++)
				for (int h=0; h<input.vList[i].targetCount.get(j); h++)
					vList[i].addTarget(input.vList[i].targets.get(j));
			for (int j=0; j<input.vList[i].sources.size(); j++)
				for (int h=0; h<input.vList[i].sourceCount.get(j); h++)
					vList[i].addSource(input.vList[i].sources.get(j));
		}
	}
	/**
	 * This constructor creates an graph copy of nodes, without a single edge.
	 * @param input FileReader
	 */
	public Graph (Graph input, int types) { //additional parameters required to manipulate group numbership
		listIndex2Id = input.listIndex2Id;
		listType2Value = input.listType2Value;
		listIndex2Type = input.listIndex2Type;
		listIndex2Value2 = input.listIndex2Value2;
		numEgs = 0;
		numVtx = input.numVtx;
		numType = input.numType;
		numTrueType = input.numTrueType;
		directed = input.directed;
		selfLoop = input.selfLoop;
		multiEdge = input.multiEdge;
		degreeCorrect = input.degreeCorrect;
		typeP = input.typeP;
		gNode = input.gNode;
		
		vtxSelfloopCount=new int[numVtx];
		vList = new Vertex[numVtx];
		for (int i=0; i<numVtx; i++) {
			vtxSelfloopCount[i]=input.vtxSelfloopCount[i];
			vList[i] = new Vertex(); //initialize the edge lists
			vList[i].index = input.vList[i].index;
			vList[i].id = input.vList[i].id;
			vList[i].type = input.vList[i].type;
			vList[i].value = input.vList[i].value;
		}
	}
	// --- Instance Methods ------------------------------------------------------ 
	/**
	 * This method removes a vertex and all its incident edges from the graph.
	 * @param index int
	 */	
	public void removeNode(int index) {
		for (int i=0; i<numVtx; i++) {
			if (i!=index) {	
				vList[i].targets.remove(Integer.valueOf(index));
				vList[i].sources.remove(Integer.valueOf(index));
			}
		}
		vList[index].targets.clear();
		vList[index].targetCount.clear();
		vList[index].outDegree = 0;
		vList[index].sources.clear();
		vList[index].sourceCount.clear();
		vList[index].inDegree = 0;
	}
	/**
	 * This method calculates log factorial. 
	 * @param t int
	 */	
	private double poisson(double base, int pow) {
		if (pow==0)
			return 1;
		double temp = 1;
		for (int i=1; i<=pow; i++)
			temp = temp*base/i;
		return temp;
	}
	/**
	 * This method draws a Poisson random variable
	 * @param t int
	 */	
	private int poissonDraw(double lambda, Random r) {
		int k = -1;
		double cdf = r.nextDouble();
		while (cdf >= 0) {
			k++;
			cdf -= poisson(lambda, k) * java.lang.Math.exp(-lambda);
		}
		return k;	
	}
		
	/**
	 * This method draws a Bernoulli random variable
	 * @param t int
	 */	
	private int bernoulliDraw(double p, Random r) {
		int k = 0;
		double temp = r.nextDouble();		
		if (temp < p)
			k++;	
		return k;	
	}
	/**
	 * This method randomize edges according the Poisson block model
	 * @param original GraphBP
	 */	
	public Graph RandomizeEdge(boolean DC, boolean newDC) {
		Graph newGraph = new Graph(this, this.getNumTrueType());
		newGraph.degreeCorrect = newDC;
		/*
		Random r = new Random(); //group assignment randomization
		for (int i=0; i<getNumNodes(); i++) {
			if (r.nextDouble() < 0.5)
				newGraph.vList[i].type = 0;
			else 
				newGraph.vList[i].type = 1;
		}
		*/
		
		Random r = new Random();
		
		for (int i=0; i<getNumNodes(); i++) {
			int degree = 1; 
			if (i<getNumNodes()/2)
				degree = 3;
			for (int j=0; j<getNumNodes(); j++) {
				int degree2 = 1;
				if (j<getNumNodes()/2)
					degree2 = 3;
				if (this.isDirected()) {
					if (j!=i) {
						int times = bernoulliDraw(typeP[newGraph.vList[i].type][newGraph.vList[j].type], r);
						if (DC)
							times = bernoulliDraw(degree*degree2*typeP[newGraph.vList[i].type][newGraph.vList[j].type], r);
						for (int h=0; h<times; h++) {
							if (newGraph.vList[i].addTarget(j))
								newGraph.numEgs++; //number of edges		;
							newGraph.vList[j].addSource(i);
						}
					}
				}
				else {
					if (j>i) {
						int times = bernoulliDraw(typeP[newGraph.vList[i].type][newGraph.vList[j].type], r);
						if (DC)
							times = bernoulliDraw(degree*degree2*typeP[newGraph.vList[i].type][newGraph.vList[j].type], r);
						for (int h=0; h<times; h++) {
							if (newGraph.vList[i].addTarget(j))
								newGraph.numEgs++; //number of edges		;
							newGraph.vList[j].addSource(i);
							newGraph.vList[j].addTarget(i);
							newGraph.vList[i].addSource(j);
						}
					}
				}
			}
		}
		return newGraph;
	}
	
	/**
	 * This method randomize edges according the vanilla Bernoulli block model
	 * @param original GraphBP
	 */	
	public Graph RandomizeEdge() {
		Graph newGraph = new Graph(this, this.getNumTrueType());
		newGraph.degreeCorrect = false;
		
		Random r = new Random(); 
		for (int i=0; i<getNumNodes(); i++) {
			for (int j=0; j<getNumNodes(); j++) {
				if (this.isDirected()) {
					if (j!=i) {
						double p = r.nextDouble();
						if (p < typeP[vList[i].type][vList[j].type]) {
							if (newGraph.vList[i].addTarget(j))
								newGraph.numEgs++; //number of edges		;
							newGraph.vList[j].addSource(i);
						}
						
					}
				}
				else {
					if (j>i) {
						double p = r.nextDouble();
						if (p < typeP[vList[i].type][vList[j].type]) {
							if (newGraph.vList[i].addTarget(j))
								newGraph.numEgs++; //number of edges		;
							newGraph.vList[j].addSource(i);
							newGraph.vList[j].addTarget(i);
							newGraph.vList[i].addSource(j);
						}
					}
				}
			}
		}
		return newGraph;
	}
	
	public void growRandEdge(double[] n, double[][] groupP) {
		Random r = new Random(); 
		for (int i=0; i<getNumNodes(); i++) {
			for (int j=0; j<getNumNodes(); j++) {
				if (this.isDirected()) {
					if (j!=i) {
						double p = r.nextDouble();
						if (p < groupP[vList[i].type][vList[j].type]) {
							if (vList[i].addTarget(j))
								numEgs++; //number of edges		;
							vList[j].addSource(i);
						}
					}
				}
				else {
					if (j>i) {
						double p = r.nextDouble();
						if (p < groupP[vList[i].type][vList[j].type]) {
							if (vList[i].addTarget(j))
								numEgs++; //number of edges		;
							vList[j].addSource(i);
							vList[j].addTarget(i);
							vList[i].addSource(j);
						}
					}
				}
			}
		}
	}
	
	/**
	 * This method consolidate the graph by removing all isolated singletons
	 * @param index int
	
	public void consolidate() {
		
	} */	
	
	/**
	 * This method returns the number of vertices.
	 * @param null
	 */	
	public int getNumNodes() { return numVtx; }
	
	/**
	 * This method returns the number of edges.
	 * @param null
	 */	
	public int getNumEdgs() { return numEgs; }
	
	/**
	 * This method checks if edges exist between given vertices.
	 * @param null
	 */	
	public boolean checkEdge(int source, int target) {
		return vList[source].targets.contains(target);
	}	
	/**
	 * This method checks if the graph is directed.
	 * @param null
	 */		
	public boolean isDirected(){ return directed; }
	/**
	 * This method checks if the graph has self-loop.
	 * @param null
	 */		
	public boolean hasSelfloop(){ return selfLoop; }	
	/**
	 * This method returns the number of true types in the graph.
	 * @param null
	 */			
	public int getNumTrueType(){ return numTrueType; }
	/**
	 * This method returns the number of types in the graph.
	 * @param null
	 */			
	public int getNumType(){ return numType; }
	/**
	 * This method returns the number of nodes in each type
	 * @param null
	 */			
	public int[] getNumNodesType() {
		int[] count = new int[numTrueType];
		for (int i=0; i<numVtx; i++)
			count[vList[i].type]++;
		return count; 
	}
	/**
	 * This method returns the number of edges between types
	 * @param null
	 */			
	public int[][] getNumEdgsType() {
		int[][] count = new int[numTrueType][numTrueType];
		for (int i=0; i<numVtx; i++)
			for (int j=0; j<numVtx; j++) if (j!=i)
				if (vList[i].targets.contains(j))
					count[vList[i].type][vList[j].type] += vList[i].targetCount.get(vList[i].targets.indexOf(j));
		return count; 	
	}
	/**
	 * This method returns the number of degree-corrected edges between types
	 * @param null
	 */			
	public double[][] getNumDegreeType() {
		double[][] count = new double[numTrueType][numTrueType];
		for (int i=0; i<numVtx; i++) {
			double iDegree = vList[i].outDegree;
			for (int j=0; j<numVtx; j++) if (j!=i) {
				double jDegree = vList[j].inDegree;
				if (vList[i].targets.contains(j))
					count[vList[i].type][vList[j].type] += (double)(vList[i].targetCount.get(vList[i].targets.indexOf(j)) / iDegree / jDegree);
			}
		}
		return count; 	
	}
	/**
	 * This method returns the lowest out degree of a node in graph
	 * @param null
	 */			
	public int getOutDegreeBase() {
		int count = Integer.MAX_VALUE;
		for (int i=0; i<numVtx; i++) {
			int temp = vList[i].outDegree;
			if (temp<count)
				count=temp;
		}
		return count; 	
	}
	/**
	 * This method returns the lowest in degree of a node in graph
	 * @param null
	 */			
	public int getInDegreeBase() {
		int count = Integer.MAX_VALUE;
		for (int i=0; i<numVtx; i++) {
			int temp = vList[i].inDegree;
			if (temp<count)
				count=temp;
		}
		return count; 	
	}
	/**
	 * This method returns the highest out degree of a node in graph
	 * @param null
	 */			
	public int getOutDegreeTop() {
		int count = 0;
		for (int i=0; i<numVtx; i++) {
			int temp = vList[i].outDegree;
			if (temp>count)
				count=temp;
		}
		return count; 	
	}
	/**
	 * This method returns the highest in degree of a node in graph
	 * @param null
	 */			
	public int getInDegreeTop() {
		int count = 0;
		for (int i=0; i<numVtx; i++) {
			int temp = vList[i].inDegree;
			if (temp>count)
				count=temp;
		}
		return count; 	
	}
}
