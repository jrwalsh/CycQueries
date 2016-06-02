package edu.iastate.cycqueries.queryScripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import edu.iastate.javacyco.Frame;
import edu.iastate.javacyco.JavacycConnection;
import edu.iastate.javacyco.PtoolsErrorException;


/**
 * Written for Liam of the CBiRC group.
 * 
 * Designed to query all genes and retrieve their GO annotations (from their product) separated by GO category.
 * Returns a tab delimited file with FrameID, GO-BiologicalProcess, GO-MolecularFunction, GO-CellularComponent, GO-Other
 */
public class LiamGOTerms {
	private static String host = "localhost";
	private static String organism = "ECOLI";
	private static int port = 4444;
	private static String geneFile = "/home/jesse/Desktop/geneList.csv";
	
	
	public static void main(String[] args) {
		try {
			Long start = System.currentTimeMillis();
			queryGOAnnotations();
			Long stop = System.currentTimeMillis();
			Long runtime = (stop - start) / 1000;
			System.out.println("Runtime is " + runtime + " seconds.");
		}
		catch(Exception e) {
			e.printStackTrace();
			System.out.println("Caught a "+e.getClass().getName()+". Shutting down...");
		}
	}
	
	static private void queryGOAnnotations() throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port);
		conn.selectOrganism(organism);
		ArrayList<Frame> allFrames = readFrameIDs(geneFile);
		
		String results = "";
		int numColumns = 5;
		results = pushTabbedLine("FrameID\tGO-BiologicalProcess\tGO-MolecularFunction\tGO-CellularComponent\tGO-Other", results, numColumns);
		
		if (allFrames == null || allFrames.isEmpty()) {
			System.err.println("No frames returned, exiting program.");
			return;
		}
		
		for (Frame gene : allFrames) {
			// For each frame, get its products
			ArrayList<Frame> products = getProductsOfGene(gene.getLocalID(), conn);

			if (products == null || products.isEmpty()) results = pushTabbedLine(gene.getLocalID(), results, numColumns);
			else {
				
				String GObp = "\t";
				String GOmf = "\t";
				String GOcc = "\t";
				String GOunk = "\t";
				
//				if (products.size() > 1) System.err.println("Here: " + gene.getLocalID());
				for (Frame product : products) {
					// For each product, get the GO-TERMS slot
					ArrayList<String> goTerms = getGOTermsOfProtein(product.getLocalID(), conn);
					
//					if (goTerms == null || goTerms.isEmpty()) results = pushTabbedLine(gene.getLocalID(), results, numColumns);
//					else {
					for (String goTerm : goTerms) {
						if (Frame.load(conn, goTerm).isGFPClass("|GO:0008150|")) { //Child of GO:0008150 is a Biological Process
							GObp += goTerm.replace("|", "") + "-" + Frame.load(conn, goTerm).getCommonName() + ",";
						} else if (Frame.load(conn, goTerm).isGFPClass("|GO:0003674|")) { //Child of GO:0003674 is a Molecular Function
							GOmf += goTerm.replace("|", "") + "-" + Frame.load(conn, goTerm).getCommonName() + ",";
						} else if (Frame.load(conn, goTerm).isGFPClass("|GO:0005575|")) { //Child of GO:0005575 is a Cellular Component
							GOcc += goTerm.replace("|", "") + "-" + Frame.load(conn, goTerm).getCommonName() + ",";
						} else {
//								System.err.println("Unknown GO-Term type for term : " + goTerm);
							GOunk += goTerm.replace("|", "") + "-" + Frame.load(conn, goTerm).getCommonName() + ",";
						}
					}
//						for (String goTerm : goTerms) {
//							// For each GO-TERM, get the CITATIONS annotation
//							ArrayList<String> citations = getCitationsOfGOTerm(goTerm, product.getLocalID(), conn);
//							
//							if (citations == null || citations.isEmpty()) results = pushTabbedLine(gene.getLocalID() + "\t" + gene.getCommonName() + "\t" + product.getLocalID() + "\t" + goTerm, results, 5);
//							else {
//								for (String citation : citations) {
//									results = pushTabbedLine(gene.getLocalID() + "\t" + gene.getCommonName() + "\t" + product.getLocalID() + "\t" + goTerm + "\t" + citation, results, 5);
//								}
//							}
//						}
//					}
				}
				if (GObp.length() > 0) GObp = GObp.substring(0,GObp.length() -1);
				if (GOmf.length() > 0) GOmf = GOmf.substring(0,GOmf.length() -1);
				if (GOcc.length() > 0) GOcc = GOcc.substring(0,GOcc.length() -1);
				if (GOunk.length() > 0) GOunk = GOunk.substring(0,GOunk.length() -1);
				results = pushTabbedLine(gene.getLocalID() + GObp + GOmf + GOcc + GOunk, results, numColumns);
			}
		}
		
		System.out.println(results);
	}
	
	// Simple Queries
	static private ArrayList<Frame> getProductsOfGene(String geneID, JavacycConnection conn) {
		ArrayList<Frame> products = new ArrayList<Frame>();
		try {
//			Gene gene = (Gene) Gene.load(conn, geneID);
//			products = gene.getSlotValues("PRODUCT");
			for (Object productID : conn.allProductsOfGene(geneID)) products.add(Frame.load(conn, (String) productID));
		} catch (PtoolsErrorException e) {
			e.printStackTrace();
		}
		return products;
	}
	
	@SuppressWarnings("unchecked")
	static private ArrayList<String> getGOTermsOfProtein(String proteinID, JavacycConnection conn) {
		ArrayList<String> goTerms = new ArrayList<String>();
		try {
			Frame protein = Frame.load(conn, proteinID);
			goTerms = protein.getSlotValues("GO-TERMS");
		} catch (PtoolsErrorException e) {
			e.printStackTrace();
		}
		return goTerms;
	}
	
//	static private ArrayList<String> getCitationsOfGOTerm(String goTerm, String proteinID, JavacycConnection conn) {
//		ArrayList<String> citations = new ArrayList<String>();
//		try {
//			citations = conn.getValueAnnots(proteinID, "GO-TERMS", goTerm, "CITATIONS");
//		} catch (PtoolsErrorException e) {
//			e.printStackTrace();
//		}
//		return citations;
//	}
	
	
	// Utilities
	static private ArrayList<Frame> readFrameIDs(String fileName) {
		// Should be csv file with column one as frameID of genes.  Ignores other columns, assumes headers and ignores them.
		JavacycConnection conn = new JavacycConnection("localhost",4444);
		conn.selectOrganism("ECOLI");
		
		ArrayList<Frame> geneFrames = new ArrayList<Frame>();
		
		File geneFile = new File(fileName);
		BufferedReader reader = null;
		
		try {
			reader = new BufferedReader(new FileReader(geneFile));
			String text = null;
			reader.readLine();
			
			while ((text = reader.readLine()) != null) {
				String frameID = text.split(",")[0];
				geneFrames.add(Frame.load(conn, frameID));
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
		finally {
			try {
			} catch (Exception e) {
				e.printStackTrace();
			}
			try {
				if (reader != null) {
					reader.close();
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		return geneFrames;
	}
	
	static private String pushTabbedLine(String newLine, String currentString, int expectedColumns) {
		String returnString = currentString;
		String addLine = newLine;
		
		int count = 0;
		for (int i = 0; i < newLine.length(); i++) {
			char c = newLine.charAt(i);
			if (c == '\t') count++;
		}
		
		int tabsToAdd = expectedColumns-count-1; // For i columns, we need i-1 tabs
		while (tabsToAdd > 0) {
			addLine += "\t";
			tabsToAdd--;
		}
		addLine += "\n";
		
		return returnString += addLine;
	}
}
