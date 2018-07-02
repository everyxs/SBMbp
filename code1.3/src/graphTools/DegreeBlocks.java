package graphTools;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Random;

/**
 * This class generates synthetic networks according to a degree-corrected block model (degree-sequence manually set)
 * The output file is in GML format, which can be directly fed into the Graph class
 * The main function has been set to produce one of the synthetic graphs used in paper
 * 
 * @author Xiaoran Yan ( xiaoran.a.yan@gmail.com )
 * @version BP_1.3
 * @time Jun, 2018
 */
public class DegreeBlocks {
	
	String path;  //output directory
	int n; //number of nodes
	int k; //number of blocks
		
	// --- Constructors ---------------------------------------------------------- 
	public DegreeBlocks(int v, int g, String dir, double assort, double disassort, int power) {
		n = v;
		k = g;
		path = dir; //output directory
		FileOutputStream output = null;
		
		String out = path.concat("/test"+n); //index the output files
		try {
			output = new FileOutputStream(out + "g5DC.gml"); //output file suffix
		}
		catch (FileNotFoundException ex) {
			ex.printStackTrace();
		}
		PrintStream print = new PrintStream(output);
		
		print.println("graph"); //start of the GML output
		print.println("[");
		
		print.println("  directed 0"); //directness of the graph
		for (int i=0; i<n; i++) {
			print.println("  node");
			print.println("  [");
			print.println("    id "+i); //node id
			print.println("    value "+ i%k); //ground truth label
			print.println("  ]");
		}
		
		//Random edge generation according to the block parameters
		Random r = new Random();
		for (int i=0; i<n; i++) {
			int degree = 3; //out-degree for half of the nodes
			if (i<n/2)
				degree = 3*power; //power times out-degree for the other half of the nodes
			for (int j=i+1; j<n; j++) {
				int degree2 = 3; //in-degree for half of the nodes
				if (j<n/2)
					degree2 = 3*power; //power times in-degree for the other half of the nodes
				for (int d=0; d<degree*degree2; d++) { //intra block connections
					if (i%k == j%k) {
						if (r.nextDouble() < assort) {
							print.println("  edge");
							print.println("  [");
							print.println("    source "+i);
							print.println("    target "+j);
							print.println("  ]");
						}	
					}
					else { //inter block connections
						if (r.nextDouble() < disassort) {
							print.println("  edge");
							print.println("  [");
							print.println("    source "+i);
							print.println("    target "+j);
							print.println("  ]");
						}
					}
				}
			}
		
		}
		print.println("  ]");
	}
	
	public static void main ( String[] args ) {
		String currentDir = System.getProperty("user.dir");
		String dir = currentDir + "/data/test";
		DegreeBlocks test1 = new DegreeBlocks(1000, 5, dir, 0.02, 0.002, 2); //set the block parameters here
	}
}