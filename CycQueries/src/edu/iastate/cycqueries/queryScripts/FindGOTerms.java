package edu.iastate.cycqueries.queryScripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.ArrayList;

import javax.swing.table.DefaultTableModel;

import edu.iastate.javacyco.Frame;
import edu.iastate.javacyco.JavacycConnection;
import edu.iastate.javacyco.PtoolsErrorException;


/**
 * Written for Taner Sen of the MaizeGDB/CornCyc group.
 * 
 * Given a list of go term annotations we expect to find in the database, validate that they really are there.
 */
public class FindGOTerms {
	private static String host = "jrwalsh.student.iastate.edu";
	private static String organism = "MAIZE";
	private static int port = 4444;
	private static String frameType = "|All-Genes|";
//	private static String outputFile = "/home/jesse/Desktop/go_output.txt";
	private static String outputFile = "C:\\Users\\Jesse\\Desktop\\output\\go_term_match_report.txt";
	private static String inputFile = "C:\\Users\\Jesse\\Desktop\\MatchingGO\\MaizeGDB_Annotations_Merged_May_2011_v3_tab_delim.csv";
	private static String delimiter = "\t";
	private static boolean hasHeader = true;
	
	/*
	 * Algorithm (simplified) for r85
	 * 1) Read list of go term data from file
	 * 2) Match go term data from file to entities in database (look for exact match to locus, then genemodel name, dropping the " in the database names, first exact match taken then stop looking), then get all (gene) frames of matches
	 * 3) Get all protein products of genes
	 * 4) Get all go terms for the proteins
	 * 5) Print information on matches or non-matches
	 */
	public static void main(String[] args) {
		try {
			Long start = System.currentTimeMillis();
			findGOTerms();
			Long stop = System.currentTimeMillis();
			Long runtime = (stop - start) / 1000;
			System.out.println("Runtime is " + runtime + " seconds.");
		}
		catch(Exception e) {
			e.printStackTrace();
			System.out.println("Caught a "+e.getClass().getName()+". Shutting down...");
		}
	}
	
	static private void findGOTerms() {
		JavacycConnection conn = new JavacycConnection(host, port);
//		JavacycConnection conn = new JavacycConnection(host, port, "me", "pass");
		conn.selectOrganism(organism);
		String results = "";
		
		ArrayList<Annotation> annotations = new ArrayList<Annotation>();
		String text = null;
		String line = "";
		
		// Read data from file
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(new File(inputFile)));
			if (hasHeader) reader.readLine(); //discard header
			while ((line = reader.readLine()) != null) {
				annotations.add(new Annotation(line));
			}
		} catch (FileNotFoundException exception) {
			exception.printStackTrace();
		} catch (Exception exception) {
			exception.printStackTrace();
		}
		finally {
			try {
			} catch (Exception exception) {
				exception.printStackTrace();
			}
			try {
				if (reader != null) {
					reader.close();
				}
			} catch (IOException exception) {
				exception.printStackTrace();
			}
		}
		
		// Process data, compare to database info
		try {
			for (Annotation annotation : annotations) {
				Frame match = null;
				ArrayList<Frame> potentialMatches = conn.search(annotation.getLocus(), frameType);
				for (Frame potentialMatch : potentialMatches) {
					ArrayList<String> frameNames = potentialMatch.getNames();
					ArrayList<String> frameNamesMod = new ArrayList<String>();
					for (String name : frameNames) frameNamesMod.add(name.replace("\"", ""));
					if (frameNamesMod.contains(annotation.getLocus())) {
						match = potentialMatch;
						break;
					}
				}
				
				if (match == null) {
					potentialMatches = conn.search(annotation.getGeneModel(), frameType);
					for (Frame potentialMatch : potentialMatches) {
						ArrayList<String> frameNames = potentialMatch.getNames();
						ArrayList<String> frameNamesMod = new ArrayList<String>();
						for (String name : frameNames) frameNamesMod.add(name.replace("\"", ""));
						if (frameNamesMod.contains(annotation.getGeneModel())) {
							match = potentialMatch;
							break;
						}
					}
				}
				
				// If we have a matching frame (gene/locus), then we should check if the go-annotation matches
				if (match == null) System.out.println("No match for " + annotation.getLocus());
				else {
					boolean matchGOTerm = false;
					boolean matchEV = false;
					Frame bestMatchGO = null; //first go match, or first go match with matching pmid if such exists
					Frame bestMatchEV = null;//first pmid and go match
					String EVMatchAnnotation = "";
					
					ArrayList<Frame> products = getProductsOfGene(match.getLocalID(), conn);
					for (Frame product : products) {
						// For each product, get the GO-TERMS slot
						ArrayList<String> goTerms = getGOTermsOfProtein(product.getLocalID(), conn);
						for (String goTerm : goTerms) {
							if (goTerm.equalsIgnoreCase(annotation.getGo())) {
								matchGOTerm = true;
								if (bestMatchGO == null) {
									bestMatchGO = product;
									EVMatchAnnotation = goTerm;
								}
								// For each GO-TERM, get the CITATIONS annotation
								ArrayList<String> citations = getCitationsOfGOTerm(goTerm, product.getLocalID(), conn);
								for (String citation : citations) {
									try {
										String evCode = citation.split(":")[1];
										if (evCode.equalsIgnoreCase(annotation.getEvCode())) {
											matchEV = true;
											if (bestMatchEV == null) {
												bestMatchEV = product;
												bestMatchGO = product;
												EVMatchAnnotation = goTerm + "\t" + citation;
											}
										}
									} catch (ArrayIndexOutOfBoundsException e) {
										//ignore
									}
								}
							}
						}
					}
					
					if (matchEV) System.out.println(annotation.getLine() + "\tGO and EV match" + "\t" + bestMatchEV.getLocalID() + "\t" + EVMatchAnnotation);
					else if (matchGOTerm) System.out.println(annotation.getLine() + "\tGO match" + "\t" + bestMatchGO.getLocalID() + "\t" + EVMatchAnnotation);
					else System.out.println(annotation.getLine() + "\tNo match found");
				}
			}
		} catch (PtoolsErrorException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}	
	
	static private ArrayList<Frame> getProductsOfGene(String geneID, JavacycConnection conn) {
		ArrayList<Frame> products = new ArrayList<Frame>();
		try {
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

	static public class Annotation {
		private String line = null;
		private String geneModel = null;
		private String locusID = null;
		private String locus = null;
		private String termName = null;
		private String go = null;
		private String pmid = null;
		private String evCode = null;
		private String timeStamp = null;
		private String source = null;
		private String comment = null;
		
		public Annotation (String line) {
			processDataLine(line);
		}
		
		// Process tab delimited data
		private void processDataLine(String line) {
			this.line = line;
			String[] splitLine = line.split(delimiter);
			try {
				geneModel = splitLine[0];
				locusID = splitLine[1];
				locus = splitLine[2];
				termName = splitLine[3];
				go = splitLine[4];
				pmid = splitLine[5];
				evCode = splitLine[6];
				timeStamp = splitLine[7];
				source = splitLine[8];
				comment = splitLine[9];
			} catch (IndexOutOfBoundsException e) {
				// ignore
			}
		}
		
		public String getLine() {
			return line;
		}
		
		public String getGeneModel() {
			return geneModel;
		}

		public String getLocusID() {
			return locusID;
		}

		public String getLocus() {
			return locus;
		}

		public String getTermName() {
			return termName;
		}

		public String getGo() {
			return go;
		}

		public String getPmid() {
			return pmid;
		}

		public String getEvCode() {
			return evCode;
		}

		public String getTimeStamp() {
			return timeStamp;
		}

		public String getSource() {
			return source;
		}

		public String getComment() {
			return comment;
		}
	}
}
