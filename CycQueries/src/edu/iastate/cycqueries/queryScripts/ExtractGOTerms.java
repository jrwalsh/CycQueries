package edu.iastate.cycqueries.queryScripts;

import java.io.File;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.ArrayList;

import edu.iastate.javacyco.Frame;
import edu.iastate.javacyco.JavacycConnection;
import edu.iastate.javacyco.PtoolsErrorException;


/**
 * Written for Taner Sen of the MaizeGDB/CornCyc group.
 * 
 * Designed to query all genes and retrieve their GO annotations (from their product) and evidence codes for the GO annotations.
 * Returns a tab delimited file with FrameID, GeneCommonName, GO-Term, EvidenceCode
 */
public class ExtractGOTerms {
	private static String host = "jrwalsh-server.student.iastate.edu";
	private static String organism = "MAIZE";
	private static int port = 4444;
	private static String frameType = "|All-Genes|";
//	private static String outputFile = "/home/jesse/Desktop/go_output.txt";
	private static String outputFile = "C:\\Users\\Jesse\\Desktop\\output\\go_output.txt";
	
	/*
	 * Algorithm (simplified) for r85
	 * 1) Get all (gene) frames of type |All-Genes|
	 * 2) Get all protein products of genes
	 * 3) Get all go terms for the proteins
	 * 4) Print gene, protein, go term, and citation information to tab delimited file
	 */
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
//		JavacycConnection conn = new JavacycConnection(host, port, "me", "pass");
		conn.selectOrganism(organism);
		String results = "";
		
		ArrayList<Frame> allFrames = getAllOfFrameType(frameType, conn);
		
		if (allFrames == null || allFrames.isEmpty()) {
			System.err.println("No frames returned, exiting program.");
			return;
		}
		
		for (Frame gene : allFrames) {
			// For each frame, get its products
			ArrayList<Frame> products = getProductsOfGene(gene.getLocalID(), conn);

			if (products == null || products.isEmpty()) results = pushTabbedLine(gene.getLocalID() + "\t" + gene.getCommonName(), results, 5);
			else {
				for (Frame product : products) {
					// For each product, get the GO-TERMS slot
					ArrayList<String> goTerms = getGOTermsOfProtein(product.getLocalID(), conn);
					
					if (goTerms == null || goTerms.isEmpty()) results = pushTabbedLine(gene.getLocalID() + "\t" + gene.getCommonName() + "\t" + product.getLocalID(), results, 5);
					else {
						for (String goTerm : goTerms) {
							// For each GO-TERM, get the CITATIONS annotation
							ArrayList<String> citations = getCitationsOfGOTerm(goTerm, product.getLocalID(), conn);
							
							if (citations == null || citations.isEmpty()) results = pushTabbedLine(gene.getLocalID() + "\t" + gene.getCommonName() + "\t" + product.getLocalID() + "\t" + goTerm, results, 5);
							else {
								for (String citation : citations) {
									citation = citation.replaceAll("\"", "");
									String decodedCitation = decodeCitationTimeStamp(citation, conn);
									results = pushTabbedLine(gene.getLocalID() + "\t" + gene.getCommonName() + "\t" + product.getLocalID() + "\t" + goTerm + "\t" + decodedCitation, results, 5);
								}
							}
						}
					}
				}
			}
		}
		
		System.out.println(results);
		printString(outputFile, results);
	}
	

	// Simple Queries
	static private ArrayList<Frame> getAllOfFrameType(String type, JavacycConnection conn) {
		ArrayList<Frame> allFrames = new ArrayList<Frame>();
		try {
			allFrames = conn.getAllGFPInstances(type);
		} catch (PtoolsErrorException e) {
			e.printStackTrace();
		}
		return allFrames;
	}
	
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
	
	@SuppressWarnings("unchecked")
	static private ArrayList<String> getCitationsOfGOTerm(String goTerm, String proteinID, JavacycConnection conn) {
		ArrayList<String> citations = new ArrayList<String>();
		try {
			citations = conn.getValueAnnots(proteinID, "GO-TERMS", goTerm, "CITATIONS");
		} catch (PtoolsErrorException e) {
			e.printStackTrace();
		}
		return citations;
	}
	
	
	// Utilities
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

	private static String decodeCitationTimeStamp(String citation, JavacycConnection conn) {
		if (citation == null || citation.equalsIgnoreCase("")) return "";
		
		String parsedString = "";
		try {
			String[] citationArray = citation.split(":");
			int timeStampIndex;
			if (citationArray.length == 4) timeStampIndex = 2;
			else if (citationArray.length == 5) timeStampIndex = 3;
			else {
				return citation;
			}
			
			String encodedTime = citationArray[timeStampIndex];
			SimpleDateFormat simpleDate = new SimpleDateFormat("MM-dd-yyyy HH-mm-ss");
			String decodedTime = simpleDate.format(conn.decodeTimeStamp(encodedTime).getTime());
			
			for (int i = 0; i < citationArray.length ; i++) {
				if (i == timeStampIndex && i == citationArray.length-1) parsedString += decodedTime;
				else if (i == timeStampIndex) parsedString += decodedTime + ":";
				else if (i == citationArray.length-1) parsedString += citationArray[i];
				else parsedString += citationArray[i] + ":";
			}
		} catch (Exception e) {
			System.err.println("Error parsing citation: " + citation);
			return citation;
		}
		return parsedString;
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
