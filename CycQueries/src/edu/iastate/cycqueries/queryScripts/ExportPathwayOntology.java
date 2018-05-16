package edu.iastate.cycqueries.queryScripts;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;

import edu.iastate.javacyco.Frame;
import edu.iastate.javacyco.JavacycConnection;
import edu.iastate.javacyco.Network;
import edu.iastate.javacyco.Pathway;
import edu.iastate.javacyco.PtoolsErrorException;

/**
 * Written for Jonathan Cohn of Syngenta
 * 
 * This was a high-level request for information about the pathway ontology in CornCyc/MaizeCyc.  He wanted to have a list of "individual pathways within the higher level 
 * categories for Zea Mays".  Here we provide 2 methods, one to get the entire pathway structure as a graph file (.gml).  Second, to sort the pathways into their top-level
 * parent classes (the ones under |Pathways|).
 * 
 * @author Jesse R. Walsh (MaizeGDB)
 */
public class ExportPathwayOntology {
	private static String host = "localhost";
	private static String organismCorn = "CORN";
	private static int port = 4444;
	
	public static void main(String[] args) {
		JavacycConnection connection = new JavacycConnection(host,port);
		connection.selectOrganism(organismCorn);
		
//		writeNetworkGML(connection);
		writePathwayTopClasses(connection);
	}
	
	public static void writeNetworkGML(JavacycConnection connection) {
		try {
			Network hierarchy = connection.getClassHierarchy(Pathway.GFPtype, true, true);
			hierarchy.writeGML(new PrintStream(new File("CornCyc_Pathway_Structure.gml")), false,false,true,false,false,false);
		} catch (PtoolsErrorException e) {
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	public static void writePathwayTopClasses(JavacycConnection connection) {
		ArrayList<String> allInstances;
		String[] classes = new String[]{
				"|Activation-Inactivation-Interconversion|", 
				"|Biosynthesis|", 
				"|Degradation|",
				"|Detoxification|",
				"|Energy-Metabolism|",
				"|Macromolecule-Modification|",
				"|Metabolic-Clusters|",
				"|Super-Pathways|",
				"|Transport-Pathways|"
				};
		String out = "TopClass\tFrameID\tCommonName\n";
		try {
			for (String thisClass : classes) {
				allInstances = connection.getClassAllInstances(thisClass);
				for (String frameName : allInstances) {
					Frame frame = Frame.load(connection, frameName);
					out += thisClass + "\t" + frame.getLocalID() + "\t" + frame.getCommonName() + "\n";
				}
			}
			printString("pathway_top_classes_CornCyc_v8.tab", out);
		} catch (PtoolsErrorException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 * Simple function to print a string to the specified file location.
	 * 
	 * @param fileName
	 * @param printString
	 */
	private static void printString(String fileName, String printString) {
		PrintStream o = null;
		try {
			o = new PrintStream(new File(fileName));
			o.println(printString);
			o.close();
		}
		catch(Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
	}
}
