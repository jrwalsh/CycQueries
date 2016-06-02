package edu.iastate.cycqueries.queryScripts;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import edu.iastate.javacyco.Compound;
import edu.iastate.javacyco.Frame;
import edu.iastate.javacyco.Gene;
import edu.iastate.javacyco.JavacycConnection;
import edu.iastate.javacyco.Network;
import edu.iastate.javacyco.Pathway;
import edu.iastate.javacyco.PtoolsErrorException;
import edu.iastate.javacyco.Reaction;

/**
 * Written for Taner Sen of the MaizeGDB/CornCyc group.
 * 
 * 
 */
public class CompareCornCycToMaizeCyc {
	private static String host = "jrwalsh.student.iastate.edu";
	private static String organismCorn = "CORN";
	private static String organismMaize = "MAIZE";
	private static int port = 4444;
	
	public static void main(String[] args) {
		try {
			Long start = System.currentTimeMillis();
			
			// Look at a specific frame in detail
			printFrame();
			
			// Just get a feel for what types of information I should look at for each frame type
//			getFrameByName(organismCorn, "|All-Genes|");
//			getFrameByName(organismCorn, "|Compounds|");
//			getFrameByName(organismCorn, "|Reactions|");
//			getFrameByName(organismCorn, "|Pathways|");
			
			
			// Collect info on specific frame times for CornCyc 5.0 and MaizeCyc 2.2
//			printGenesTab(organismCorn, "geneCounts_CornCyc5.tab");
//			printGenesTab(organismMaize, "geneCounts_MaizeCyc2.tab");
//			
//			printCompoundsTab(organismCorn, "compoundCounts_CornCyc5.tab");
//			printCompoundsTab(organismMaize, "compoundCounts_MaizeCyc2.tab");
//			
//			printReactionsTab(organismCorn, "reactionCounts_CornCyc5.tab");
//			printReactionsTab(organismMaize, "reactionCounts_MaizeCyc2.tab");
//			
//			printPathwaysTab(organismCorn, "pathwayCounts_CornCyc5.tab");
//			printPathwaysTab(organismMaize, "pathwayCounts_MaizeCyc2.tab");
			
			
			// Collect info on specific frame times for CornCyc 4.0.1
//			printGenesTab(organismCorn, "geneCounts_CornCyc4.tab");
//			printCompoundsTab(organismCorn, "compoundCounts_CornCyc4.tab");
//			printReactionsTab(organismCorn, "reactionCounts_CornCyc4.tab");
//			printPathwaysTab(organismCorn, "pathwayCounts_CornCyc4.tab");
			
			
			
			// This gets the counts to compare against Jackies counts.
//			getAllCounts(organismMaize);
//			getAllCounts(organismCorn);
			
			Long stop = System.currentTimeMillis();
			Long runtime = (stop - start) / 1000;
			System.out.println("Runtime is " + runtime + " seconds.");
		}
		catch(Exception e) {
			e.printStackTrace();
			System.out.println("Caught a "+e.getClass().getName()+". Shutting down...");
		}
	}
	
	static private void getAllCounts(String organism) throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port);
		conn.selectOrganism(organism);
		
		getGeneCounts(conn, organism);
		getCompoundCounts(conn, organism);
		getReactionCounts(conn, organism);
		getEnzymaticReactionCounts(conn, organism);
		getPathwayCounts(conn, organism);
//		getGOAnnotationCounts(conn, organism);
	}
	
	static private void getGeneCounts(JavacycConnection conn, String organism) throws PtoolsErrorException {
		conn.selectOrganism(organism);
		
		Network geneHierarchy = conn.getClassHierarchy("|Genes|", true, true);
		Set<Frame> geneNodes = geneHierarchy.getNodes();
		
		HashSet<String> geneClasses = new HashSet<String>();
		HashSet<String> geneInstances = new HashSet<String>();
		HashSet<String> geneUniqueClasses = new HashSet<String>();
		HashSet<String> geneUniqueInstances = new HashSet<String>();
		for (Frame gene : geneNodes) {
			if (gene.isClassFrame()) {
				geneClasses.add(gene.getLocalID());
				geneUniqueClasses.add(gene.getCommonName().replaceAll("_P..$|_T..$", ""));
			}
			else {
				geneInstances.add(gene.getLocalID());
				geneUniqueInstances.add(gene.getCommonName().replaceAll("_P..$|_T..$", ""));
			}
		}
		System.out.println("Gene Classes: " + geneClasses.size());
		System.out.println("Gene Instances: " + geneInstances.size());
		System.out.println("Gene Unique Classes without suffix: " + geneClasses.size());
		System.out.println("Gene Unique Instances without suffix: " + geneUniqueInstances.size());
	}
	
	
	static private void getCompoundCounts(JavacycConnection conn, String organism) throws PtoolsErrorException {
		conn.selectOrganism(organism);
		
		Network compoundHierarchy = conn.getClassHierarchy("|Compounds|", true, true);
		Set<Frame> compoundNodes = compoundHierarchy.getNodes();
		
		HashSet<String> compoundClasses = new HashSet<String>();
		HashSet<String> compoundInstances = new HashSet<String>();
		for (Frame compound : compoundNodes) {
			if (compound.isClassFrame()) compoundClasses.add(compound.getLocalID());
			else compoundInstances.add(compound.getLocalID());
		}
		System.out.println("Compound Classes: " + compoundClasses.size());
		System.out.println("Compound Instances: " + compoundInstances.size());
	}
	
	
	static private void getReactionCounts(JavacycConnection conn, String organism) throws PtoolsErrorException {
		conn.selectOrganism(organism);
		
		Network reactionsHierarchy = conn.getClassHierarchy("|Reactions|", true, true);
		Set<Frame> reactionNodes = reactionsHierarchy.getNodes();

		HashSet<String> reactionClasses = new HashSet<String>();
		HashSet<String> reactionInstances = new HashSet<String>();
		for (Frame reaction : reactionNodes) {
			if (reaction.isClassFrame()) reactionClasses.add(reaction.getLocalID());
			else reactionInstances.add(reaction.getLocalID());
		}
		System.out.println("Reaction Classes: " + reactionClasses.size());
		System.out.println("Reaction Instances: " + reactionInstances.size());
	}
	
	
	static private void getEnzymaticReactionCounts(JavacycConnection conn, String organism) throws PtoolsErrorException {
		conn.selectOrganism(organism);
		
		Network enzyRxnsHierarchy = conn.getClassHierarchy("|Enzymatic-Reactions|", true, true);
		Set<Frame> enzyRxnNodes = enzyRxnsHierarchy.getNodes();
		
		HashSet<String> enzyRxnClasses = new HashSet<String>();
		HashSet<String> enzyRxnInstances = new HashSet<String>();
		for (Frame enzyRxn : enzyRxnNodes) {
			if (enzyRxn.isClassFrame()) enzyRxnClasses.add(enzyRxn.getLocalID());
			else enzyRxnInstances.add(enzyRxn.getLocalID());
		}
		System.out.println("Enzymatic-Reactions Classes: " + enzyRxnClasses.size());
		System.out.println("Enzymatic-Reactions Instances: " + enzyRxnInstances.size());
	}
	
	
	static private void getPathwayCounts(JavacycConnection conn, String organism) throws PtoolsErrorException {
		conn.selectOrganism(organism);
		
		Network pathwaysHierarchy = conn.getClassHierarchy("|Pathways|", true, true);
		Set<Frame> pathwayNodes = pathwaysHierarchy.getNodes();
		
		HashSet<String> pathwayClasses = new HashSet<String>();
		HashSet<String> superpathways = new HashSet<String>();
		HashSet<String> pathwayInstances = new HashSet<String>();
		for (Frame pathway : pathwayNodes) {
			try {
				if (pathway.isClassFrame()) pathwayClasses.add(pathway.getLocalID());
//				else if (pathway.getSlotValues("Sub-Pathways").size() > 0) superpathways.add(pathway.getLocalID()); //not all super pathways have items in the sub-pathways slot
				else if (pathway.isGFPClass("|Super-Pathways|")) superpathways.add(pathway.getLocalID());
				else pathwayInstances.add(pathway.getLocalID());
			} catch (Exception e) {
				System.out.println("Problem with pathway : " + pathway.getLocalID() + " : " + pathway.isClassFrame());
			}
		}
		System.out.println("Pathway Classes: " + pathwayClasses.size());
		System.out.println("Pathway Instances: " + pathwayInstances.size());
		System.out.println("Superpathways: " + superpathways.size());
	}
	
	static private void getGOAnnotationCounts(JavacycConnection conn, String organism) throws PtoolsErrorException {
		conn.selectOrganism(organism);
		
		
		//TODO what to count, exactly?  assignments? citations? unique or not? etc...
		Network pathwaysHierarchy = conn.getClassHierarchy("|Polypeptides|", true, true);
		Set<Frame> proteinNodes = pathwaysHierarchy.getNodes();
		
		
//		HashSet<String> gotermAssignments = new HashSet<String>();
//		HashSet<String> gotermCitations = new HashSet<String>();
		for (Frame protein : proteinNodes) {
			try {
				System.out.println(protein.getLocalID());
				for (Object goTerm : protein.getSlotValues("GO-TERMS")) {
					System.out.println("\t"+goTerm.toString());
					for(Object citation : protein.getAnnotations("GO-TERMS", goTerm.toString(), "CITATIONS")) {
						System.out.println("\t\t"+citation.toString());
					}
				}
//				if (protein.isClassFrame()) proteinClasses.add(protein.getLocalID());
//				else proteinInstances.add(protein.getLocalID());
			} catch (Exception e) {
				System.out.println("Problem with protein : " + protein.getLocalID() + " : " + protein.isClassFrame());
			}
		}
//		System.out.println("Protein Classes: " + proteinClasses.size());
//		System.out.println("Protein Instances: " + proteinInstances.size());
	}
	
	
	
	static private void printGenesTab(String organism, String fileName) throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port);
		conn.selectOrganism(organism);
		
		Network hierarchy = conn.getClassHierarchy("|All-Genes|", true, true);
		Set<Frame> nodes = hierarchy.getNodes();
		
		String printString = "";
		printString += nodes.size() + "\n";
		printString += "FrameID\tCommonName\tisClass?" + "\n";
		for (Frame node : nodes) {
			printString += node.getLocalID() + "\t" + node.getCommonName() + "\t" + node.isClassFrame() + "\n";
		}
		
		printString(fileName, printString);
	}
	
	static private void printCompoundsTab(String organism, String fileName) throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port);
		conn.selectOrganism(organism);
		
		Network hierarchy = conn.getClassHierarchy("|Compounds|", true, true);
		Set<Frame> nodes = hierarchy.getNodes();
		
		String printString = "";
		printString += nodes.size() + "\n";
		printString += "FrameID\tCommonName\tInChI\tisClass?" + "\n";
		for (Frame node : nodes) {
			printString += node.getLocalID() + "\t" + node.getCommonName() + "\t" + node.getSlotValue("INCHI") + "\t" + node.isClassFrame() + "\n";
		}
		
		printString(fileName, printString);
	}
	
	static private void printReactionsTab(String organism, String fileName) throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port);
		conn.selectOrganism(organism);
		
		Network hierarchy = conn.getClassHierarchy("|Reactions|", true, true);
		Set<Frame> nodes = hierarchy.getNodes();
		
		String printString = "";
		printString += nodes.size() + "\n";
		printString += "FrameID\tCommonName\tEC\tDirection\tNumEnzyRxns\tisClass?\tReactants\tProducts" + "\n";
		for (Frame node : nodes) {
			
			String reactants = "";
			for (Object reactant : node.getSlotValues("Left")) {
				reactants += reactant + ", ";
			}
			if (reactants.length() > 0) reactants = reactants.substring(0, reactants.length()-2);
			
			String products = "";
			for (Object product : node.getSlotValues("right")) {
				products += product + ", ";
			}
			if (products.length() > 0) products = products.substring(0, products.length()-2);
			
			printString += node.getLocalID() + "\t" + node.getCommonName() + "\t" + node.getSlotValue("EC-NUMBER") + "\t" + node.getSlotValue("REACTION-DIRECTION") + "\t" + node.getSlotValues("ENZYMATIC-REACTION").size() + "\t" + node.isClassFrame() + "\t" + reactants + "\t" + products + "\n";
		}
		
		printString(fileName, printString);
	}
	
	static private void printPathwaysTab(String organism, String fileName) throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port);
		conn.selectOrganism(organism);
		
		Network hierarchy = conn.getClassHierarchy("|Pathways|", true, true);
		Set<Frame> nodes = hierarchy.getNodes();
		
		String printString = "";
		printString += nodes.size() + "\n";
		printString += "FrameID\tCommonName\tVariants?\tisSuperPathway\tisClass?\tReaction-List" + "\n";
		for (Frame node : nodes) {
			String reactionList = "";
			for (Object reaction : node.getSlotValues("Reaction-List")) {
				reactionList += reaction + ", ";
			}
			if (reactionList.length() > 0) reactionList = reactionList.substring(0, reactionList.length()-2);
			
			printString += node.getLocalID() + "\t" + node.getCommonName() + "\t" + node.getSlotValue("Variants?") + "\t" + node.isGFPClass("|Super-Pathways|") + "\t" + node.isClassFrame() + "\t" + reactionList + "\n";
		}
		
		printString(fileName, printString);
	}
	
	static private void getFrameByName(String organism, String type) throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port);
		conn.selectOrganism(organism);
		
		Network hierarchy = conn.getClassHierarchy(type, true, true);
		Set<Frame> nodes = hierarchy.getNodes();
		
		for (Frame node : nodes) {
			node.print();
			printDirectSubs(conn, node);
		}
	}
	
	static private void printDirectSubs(JavacycConnection conn, Frame frame) {
		try {
			for (Object sub : conn.getClassDirectSubs(frame.getLocalID())) {
				Frame.load(conn, sub.toString()).print();
				printDirectSubs(conn, Frame.load(conn, sub.toString()));
			}
		} catch (PtoolsErrorException e) {
			System.out.println("Didin't like : " + frame.getLocalID());
		}
	}
	
	static private void getWholeStructure() throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port);
		conn.selectOrganism(organismCorn);
		
		Network hierarchy = conn.getClassHierarchy("OCELOT-GFP::THINGS", false, true);
		hierarchy.printStructureTab();
	}
	
	static private void printFrame() throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port);
		conn.selectOrganism(organismCorn);
		
		Frame.load(conn, "|Genes|").print();
		
//		ArrayList<String> genes = conn.getClassAllInstances("|Genes|");
//		for (String gene : genes) {
//			System.out.println(Frame.load(conn, gene).getCommonName());
//		}
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
