package graphTools;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;

/**
 * This class prune out all nodes and edges except the giant component.
 * Both the input and output file should be in the simple GML format, which can be directly fed into the Graph class.
 * 
 * @author Xiaoran Yan ( everyxt@gmail.com )
 * @version BP_1.3
 * @time Nov, 2013
 */
public class PruneGML {
	
	public Graph gOriginal; //the original graph
	int[] compList; //list for updating components
	ArrayList<Integer> compCount= new ArrayList<Integer>();
	
	// --- Constructors ---------------------------------------------------------- 
	public PruneGML(String dir) {
		FileReader input = null;
		try {
			input = new FileReader(dir + "news2f.gml");
		}
		catch (FileNotFoundException ex) {
			ex.printStackTrace();
		}
		gOriginal = new Graph(input, true);
		compList = new int[gOriginal.getNumNodes()];
		for (int i=0; i<compList.length; i++)
			compList[i] = -1; //unmarked
	}
	// --- Instance Methods ------------------------------------------------------ 
	/**
	 * Flag the component membership of a node and all its descendants recursively
	 * @param root int
	 * @param flag int
	 */
	public void flagComp(int root, int flag) {
		compList[root] = flag; //flag the root
		compCount.set(flag, compCount.get(flag)+1); //increase the component size counter
		for (int i=0; i<gOriginal.vList[root].targets.size(); i++) {
			int next = gOriginal.vList[root].targets.get(i);
			if (compList[next] == -1)
				flagComp(next, flag); //flag the targets recursively
		}
		for (int i=0; i<gOriginal.vList[root].sources.size(); i++) {
			int next = gOriginal.vList[root].sources.get(i);
			if (compList[next] == -1)
				flagComp(next, flag); //flag the sources recursively
		}
	}
	
	/**
	 * Generate the giant component in simple GML format
	 * @param dir String
	 */
	public void generate (String dir){
		//Create the output files
		FileOutputStream output = null;
		try {
			output = new FileOutputStream(dir + "news2fp.gml");
		}
		catch (FileNotFoundException ex) {
			ex.printStackTrace();
		}
		PrintStream print = new PrintStream(output);
		int direct;
		if (gOriginal.isDirected())
			direct = 1;
		else direct = 0;
		print.println("directed " + direct); //copy the directness of the original graph
		

		int comp = 0;
		for (int i=0; i<compList.length; i++) if (compList[i] == -1){
			compCount.add(0);  // new component found
			flagComp(i, comp); // flag the component
			comp++; // increase the flag for the next component
		}
		
		comp = 0;
		int max = 0;
		for (int i=0; i<compCount.size(); i++) if (compCount.get(i) > max) {
			comp = i;
			max = compCount.get(i); //get the giant component
		}
			
		for (int i=0; i<gOriginal.vList.length; i++) {
			if (compList[i] == comp) { //only print the nodes in the giant component
				print.println("id " + gOriginal.vList[i].id);
				print.println("value " + gOriginal.vList[i].value);
			}
		}
		print.println("edge"); //only print the edges in the giant component
		for (int i=0; i<gOriginal.vList.length; i++) if (compList[i] == comp) {
			for (int j=0; j<gOriginal.vList[i].targets.size(); j++) if (compList[gOriginal.vList[i].targets.get(j)] == comp){
				for (int k=0; k<gOriginal.vList[i].targetCount.get(j); k++) {
					print.println("source " + gOriginal.vList[i].id);
					print.println("target " + gOriginal.listIndex2Id.get(gOriginal.vList[i].targets.get(j)));
				}
			}
		}
		System.out.println("Gml created.");
	}
	
	public static void main ( String[] args ) {
		String dir = System.getProperty("user.dir");
		dir += "/data/brown/";
		PruneGML test1 = new PruneGML(dir);
		test1.generate(dir);
	}
}