package edu.iastate.cycqueries.queryScripts;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import java.util.TreeSet;

import edu.iastate.javacyco.Frame;
import edu.iastate.javacyco.Gene;
import edu.iastate.javacyco.JavacycConnection;
import edu.iastate.javacyco.Network;
import edu.iastate.javacyco.Protein;
import edu.iastate.javacyco.PtoolsErrorException;

/**
 * Written for Taner Sen of the MaizeGDB/CornCyc group.
 * 
 * Designed to check several assumptions about the data in MaizeCyc and CornCyc, how it is stored, and how well identifiers match between both databases.
 */
public class GeneProducts {
	private static String host = "jrwalsh.student.iastate.edu";
	private static String organism = "ECOLI";
	private static int port = 4444;
	
	public static void main(String[] args) {
		try {
			Long start = System.currentTimeMillis();
			printFrame("CPLX0-7656");
//			run();
			
			Long stop = System.currentTimeMillis();
			Long runtime = (stop - start) / 1000;
			System.out.println("Runtime is " + runtime + " seconds.");
		}
		catch(Exception e) {
			e.printStackTrace();
			System.out.println("Caught a "+e.getClass().getName()+". Shutting down...");
		}
	}
	
	private static void printFrame(String frameID) throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection("jrwalsh.student.iastate.edu", 4444);
		conn.selectOrganism("ECOLI");
		
		Frame frame = Frame.load(conn, frameID);
		frame.print();
		
		for (Object s : conn.reactionsOfEnzyme("CPLX0-7656")) {
			System.out.println(s.toString());
		}
	}
	
	// Gene products are |Polypeptides| or |RNAs| only
	static private void run() throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port);
		conn.selectOrganism(organism);
		
		Network hierarchy = conn.getClassHierarchy("|Genes|", true, true);
		Set<Frame> genes = hierarchy.getNodes();
		
		for (Frame gene : genes) {
			if (!gene.isClassFrame()) {
				try {
					ArrayList<Protein> products = ((Gene)gene).getProducts();
					for (Protein product : products) {
						System.out.println(gene.getLocalID() + "\t" + product.getLocalID());
					}
				} catch (Exception e) {
					System.err.println(gene.getLocalID());
				}
			}
		}
	}
	
}
