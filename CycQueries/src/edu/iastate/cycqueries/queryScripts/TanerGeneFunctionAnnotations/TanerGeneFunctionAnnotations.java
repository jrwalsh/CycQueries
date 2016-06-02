package edu.iastate.cycqueries.queryScripts.TanerGeneFunctionAnnotations;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import edu.iastate.javacyco.Frame;
import edu.iastate.javacyco.Gene;
import edu.iastate.javacyco.JavacycConnection;
import edu.iastate.javacyco.Protein;
import edu.iastate.javacyco.PtoolsErrorException;
import edu.iastate.javacyco.Reaction;

/**
 * Taner is looking for the citations linking an enzyme to an EC number given a starting gene model name.  In corncyc and maizecyc, we search for
 * genes which match the given term, match it to enzymes, match the enzymes to enzymatic reactions (in the catalyzes slot), and pull out the
 * citation associated with those enzymatic reactions, and match them to the reaction so we can pull the EC number.  This is to be done on the
 * publically available base versions of each database, not including any schema upgrades or update propagation steps.
 * -- (i.e. corncyc 4.0 in 16.5 and maizecyc 2.2 in 17.0)
 * 
 * @throws PtoolsErrorException
 * @author Jesse Walsh 5/24/2016
 */
public class TanerGeneFunctionAnnotations {
	private static String hostCorn = "jrwalsh.student.iastate.edu";
	private static String hostMaize = "jrwalsh-server.student.iastate.edu";
	private static String organismCorn = "CORN";
	private static String organismMaize = "MAIZE";
	private static int port = 4444;
	private static String frameType = "|All-Genes|";
	private static String inFileCorn = "C:\\Users\\Jesse\\workspace\\CycQueries\\bin\\edu\\iastate\\cycqueries\\queryScripts\\TanerGeneFunctionAnnotations\\corncyc-match-TP.txt";
	private static String inFileMaize = "C:\\Users\\Jesse\\workspace\\CycQueries\\bin\\edu\\iastate\\cycqueries\\queryScripts\\TanerGeneFunctionAnnotations\\maizecyc-match-TP.txt";

	public static void main(String[] args) {
		try {
			Long start = System.currentTimeMillis();
			
			tanerEVCodesCornCyc();
			tanerEVCodesMaizeCyc();
			
			Long stop = System.currentTimeMillis();
			Long runtime = (stop - start) / 1000;
			System.out.println("Runtime is " + runtime + " seconds.");
		}
		catch(Exception e) {
			e.printStackTrace();
			System.out.println("Caught a "+e.getClass().getName()+". Shutting down...");
		}
		
	}
	
	/**
	 * Take in a list of gene model name to EC pairs (with additional columns of info) and
	 * 1) search gene model name to resolve into a gene
	 * 2) get products of gene
	 * 3) get catalyzes of product (i.e. enzymatic reactions)
	 * 4) print citations for the enzymatic reactions
	 * 5) get reaction of enzymatic reaction
	 * 6) print EC number of reaction
	 * 
	 * @throws PtoolsErrorException
	 */
	private static void tanerEVCodesCornCyc() throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(hostCorn, port);
		System.out.println("Check EV codes for Taner: CORN");
		System.out.println("Gene_Model\tCornCyc_EC\tUniProt_EC\tMatchedGeneID\tMatchedGeneName\tReactionID\tEC-Number\tCitations");
		conn.selectOrganism(organismCorn);
		
		BufferedReader reader = null;
		String line = null;
		try {
			reader = new BufferedReader(new FileReader(inFileCorn));
			reader.readLine();
			while ((line = reader.readLine()) != null) {
				String[] data = line.split("\t");
				String item = data[0];
				
				System.out.print(data[0] + "\t" + data[1] + "\t" + data[2]); //original info
				
				ArrayList<Frame> list = conn.search(item, frameType); // Search for gene
				if (list.size() != 1) {
					System.out.println();
					System.err.println("List size err " + list.size() + " for " + item);
				} else {
					Gene gene = (Gene)list.get(0);
					System.out.print("\t" + gene.getLocalID() + "\t" + gene.getCommonName()); //matched object info
					
					ArrayList<Protein> proteins = gene.getProducts();
					if (gene.getProducts().size() != 1) {
						System.out.println();
						System.err.println("List size err " + proteins.size() + " for " + item);
					} else {
						Protein protein = proteins.get(0);
						boolean firstCatalysis = true;
						if (protein.getCatalysis().size() == 0) System.out.println();
						for (Frame catalysis : protein.getCatalysis()) {
							Reaction reaction = (Reaction) Reaction.load(conn, catalysis.getSlotValue("REACTION"));
							if (firstCatalysis) {
								System.out.print("\t" + reaction.getLocalID() + "\t" + reaction.getEC()); //EC info
								firstCatalysis = false;
							} else {
								System.out.print("\t\t\t\t\t" + reaction.getLocalID() + "\t" + reaction.getEC()); //EC info
							}
							
							boolean firstCitation = true;
							for (Object citation : catalysis.getSlotValues("CITATIONS")) {
								if (firstCitation) {
									System.out.println("\t" + citation.toString()); //Citation info
									firstCitation = false;
								} else {
									System.out.println("\t\t\t\t\t\t\t" + citation.toString()); //Citation info
								}
							}
						}
					}
				}
			}
		} catch (FileNotFoundException exception) {
			exception.printStackTrace();
		} catch (IOException exception) {
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
	}
	
	/**
	 * Take in a list of gene model name to EC pairs (with additional columns of info) and
	 * 1) search gene model name to resolve into a gene
	 * 2) get products of gene
	 * 3) get catalyzes of product (i.e. enzymatic reactions)
	 * 4) print citations for the enzymatic reactions
	 * 5) get reaction of enzymatic reaction
	 * 6) print EC number of reaction
	 * 
	 * @throws PtoolsErrorException
	 */
	private static void tanerEVCodesMaizeCyc() throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(hostMaize, port);
		System.out.println("Check EV codes for Taner: MAIZE");
		System.out.println("Gene_Model\tCornCyc_EC\tUniProt_EC\tMatchedGeneID\tMatchedGeneName\tReactionID\tEC-Number\tCitations");
		conn.selectOrganism(organismMaize);
		
		BufferedReader reader = null;
		String line = null;
		try {
			reader = new BufferedReader(new FileReader(inFileMaize));
			reader.readLine();
			while ((line = reader.readLine()) != null) {
				String[] data = line.split("\t");
				String item = data[0];
				
				System.out.print(data[0] + "\t" + data[1] + "\t" + data[2]); //original info
				
				ArrayList<Frame> list = conn.search(item, frameType); // Search for gene
				if (list.size() != 1) {
					System.out.println();
					System.err.println("List size err " + list.size() + " for " + item);
				} else {
					Gene gene = (Gene)list.get(0);
					System.out.print("\t" + gene.getLocalID() + "\t" + gene.getCommonName()); //matched object info
					
					ArrayList<Protein> proteins = gene.getProducts();
					if (gene.getProducts().size() != 1) {
						System.out.println();
						System.err.println("List size err " + proteins.size() + " for " + item);
					} else {
						Protein protein = proteins.get(0);
						boolean firstCatalysis = true;
						if (protein.getCatalysis().size() == 0) System.out.println();
						for (Frame catalysis : protein.getCatalysis()) {
							Reaction reaction = (Reaction) Reaction.load(conn, catalysis.getSlotValue("REACTION"));
							if (reaction == null) continue;
							if (firstCatalysis) {
								System.out.print("\t" + reaction.getLocalID() + "\t" + reaction.getEC()); //EC info
								firstCatalysis = false;
							} else {
								System.out.print("\t\t\t\t\t" + reaction.getLocalID() + "\t" + reaction.getEC()); //EC info
							}
							
							boolean firstCitation = true;
							if (catalysis.getSlotValues("CITATIONS").size() == 0) System.out.println();
							for (Object citation : catalysis.getSlotValues("CITATIONS")) {
								if (firstCitation) {
									System.out.println("\t" + citation.toString()); //Citation info
									firstCitation = false;
								} else {
									System.out.println("\t\t\t\t\t\t\t" + citation.toString()); //Citation info
								}
							}
						}
					}
				}
			}
		} catch (FileNotFoundException exception) {
			exception.printStackTrace();
		} catch (IOException exception) {
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
	}
}
