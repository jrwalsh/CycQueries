package edu.iastate.cycqueries.queryScripts;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeSet;

import edu.iastate.javacyco.Frame;
import edu.iastate.javacyco.Gene;
import edu.iastate.javacyco.JavacycConnection;
import edu.iastate.javacyco.Protein;
import edu.iastate.javacyco.PtoolsErrorException;

/**
 * Written for Taner Sen of the MaizeGDB/CornCyc group.
 * 
 * Designed to check several assumptions about the data in MaizeCyc and CornCyc, how it is stored, and how well identifiers match between both databases.
 */
public class MaizeCycExploration {
	private static String host = "localhost";
	private static String organismCorn = "CORN";
	private static String organismMaize = "MAIZE";
	private static String geneCategory = "|All-Genes|"; //|All-Genes| or |Genes|
	private static String proteinCategory = "|Proteins|";
	private static int port = 4444;
	
	public static void main(String[] args) {
		try {
			Long start = System.currentTimeMillis();
			
			genesToProduct(organismCorn);
			autoMatchNames(organismCorn, "/home/jesse/Desktop/maizegdb_oboterms.sample.txt", 0);
			genesToProduct(organismMaize);
			autoMatchNames(organismMaize, "/home/jesse/Desktop/maizegdb_oboterms.sample.txt", 0);
			geneNameUniqueness(organismCorn);
			geneNameUniqueness(organismMaize);
			proteinNameUniqueness(organismCorn);
			proteinNameUniqueness(organismMaize);
			compareGeneFrameExistance(organismMaize, organismCorn);
			compareGeneFrameExistanceFuzzy(organismMaize, organismCorn);
			proteinNames(organismCorn);
			
			Long stop = System.currentTimeMillis();
			Long runtime = (stop - start) / 1000;
			System.out.println("Runtime is " + runtime + " seconds.");
		}
		catch(Exception e) {
			e.printStackTrace();
			System.out.println("Caught a "+e.getClass().getName()+". Shutting down...");
		}
	}
	
	// Gene products are |Polypeptides| or |RNAs| only
	static private void genesToProduct(String organism) throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port, "you", "nopass");
		conn.selectOrganism(organism);
		ArrayList<String> genesWithNoProduct = new ArrayList<String>();
		ArrayList<String> genesWithOneProduct = new ArrayList<String>();
		ArrayList<String> genesWithManyProduct = new ArrayList<String>();
		
		for (Frame f : conn.getAllGFPInstances(geneCategory)) {
			Gene gene = (Gene)f;
			if (gene.getProducts().size() == 0) genesWithNoProduct.add(gene.getLocalID());
			else if (gene.getProducts().size() == 1) genesWithOneProduct.add(gene.getLocalID());
			else genesWithManyProduct.add(gene.getLocalID());
		}
		
		System.out.println("Genes with no product: " + genesWithNoProduct.size());
		for (String gene : genesWithNoProduct) System.out.println(gene);
		System.out.println("\n\n");
		
		System.out.println("Genes with one product: " + genesWithOneProduct.size());
		for (String gene : genesWithOneProduct) System.out.println(gene);
		System.out.println("\n\n");
				
		System.out.println("Genes with many products: " + genesWithManyProduct.size());
		for (String gene : genesWithManyProduct) System.out.println(gene);
		System.out.println("\n\n");
	}
	
	static private void autoMatchNames(String organism, String fileName, int column) throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port, "you", "nopass");
		conn.selectOrganism(organism);
		
		ArrayList<String> noMatches = new ArrayList<String>();
		ArrayList<String> nonGeneMatches = new ArrayList<String>();
		HashMap<String, String> singleMatches = new HashMap<String, String>();
		HashMap<String, ArrayList<String>> multipleMatches = new HashMap<String, ArrayList<String>>();
		
		for (String name : readFile(fileName, column)) {
			ArrayList<Frame> searchHitsGenes = conn.search(name, geneCategory);
			if (searchHitsGenes.isEmpty()) {
				ArrayList<Frame> searchHitsProteins = conn.search(name, proteinCategory);
				if (searchHitsProteins.isEmpty()) noMatches.add(name);
				else nonGeneMatches.add(name);
			}
			else if (searchHitsGenes.size() == 1) singleMatches.put(name, searchHitsGenes.get(0).getLocalID());
			else {
				ArrayList<String> matches = new ArrayList<String>();
				for (Frame hit : searchHitsGenes) matches.add(hit.getLocalID());
				multipleMatches.put(name, matches);
			}
		}
		
		System.out.println("Names with no matchs: " + noMatches.size());
		for (String name : noMatches) System.out.println(name);
		System.out.println("\n\n");
				
		System.out.println("Names with non gene matchs: " + nonGeneMatches.size());
		for (String name : nonGeneMatches) System.out.println(name);
		System.out.println("\n\n");
				
		System.out.println("Names with single gene matchs: " + singleMatches.keySet().size());
		for (String name : singleMatches.keySet()) System.out.println(name + "\t" + singleMatches.get(name));
		System.out.println("\n\n");
				
		System.out.println("Names with multiple matchs: " + multipleMatches.keySet().size());
		for (String name : multipleMatches.keySet()) {
			System.out.print(name + "\t");
			for (String match : multipleMatches.get(name)) System.out.print(match + ",");
			System.out.println();
		}
		System.out.println("\n\n");
	}
	
	@SuppressWarnings("unchecked")
	static private void compareGeneFrameExistance(String maize, String corn) throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port, "you", "nopass");
		TreeSet<String> cornGeneIDs = new TreeSet<String>();
		TreeSet<String> maizeGeneIDs = new TreeSet<String>();
		
		conn.selectOrganism(corn);
		for (Frame f : conn.getAllGFPInstances(geneCategory)) {
			Gene gene = (Gene)f;
//			cornGeneIDs.add(gene.getLocalID());
//			cornGeneIDs.add(gene.getCommonName());
			cornGeneIDs.add(gene.getCommonName().replace("_P0.", "").replace("_T0.", ""));
		}
		
		conn.selectOrganism(maize);
		for (Frame f : conn.getAllGFPInstances(geneCategory)) {
			Gene gene = (Gene)f;
//			maizeGeneIDs.add(gene.getLocalID());
//			maizeGeneIDs.add(gene.getCommonName());
			maizeGeneIDs.add(gene.getCommonName().replace("_P0.", "").replace("_T0.", ""));
		}
		
		TreeSet<String> cornUniqueGenes = (TreeSet<String>) cornGeneIDs.clone();
		cornUniqueGenes.removeAll(maizeGeneIDs);
		
		TreeSet<String> maizeUniqueGenes = (TreeSet<String>) maizeGeneIDs.clone();
		maizeUniqueGenes.removeAll(cornGeneIDs);
		
		TreeSet<String> commonGenes = (TreeSet<String>) cornGeneIDs.clone();
		commonGenes.retainAll(maizeGeneIDs);
		
		System.out.println("Genes that match: " + commonGenes.size());
		for (String name : commonGenes) {
			System.out.println(name);
		}
		System.out.println("\n\n");
		
		System.out.println("CornCyc Unique Genes: " + cornUniqueGenes.size());
		for (String name : cornUniqueGenes) {
			System.out.println(name);
		}
		System.out.println("\n\n");
		
		System.out.println("MaizeCyc Unique Genes: " + maizeUniqueGenes.size());
		for (String name : maizeUniqueGenes) {
			System.out.println(name);
		}
		System.out.println("\n\n");
	}
	
	static private void compareGeneFrameExistanceFuzzy(String maize, String corn) throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port, "you", "nopass");
		ArrayList<String> cornGeneIDs = new ArrayList<String>();
		ArrayList<String> maizeGeneIDs = new ArrayList<String>();
		
		conn.selectOrganism(corn);
		for (Frame f : conn.getAllGFPInstances(geneCategory)) {
			Gene gene = (Gene)f;
//			cornGeneIDs.add(gene.getLocalID());
//			cornGeneIDs.add(gene.getCommonName());
			cornGeneIDs.add(gene.getCommonName().replaceAll("_P0.", "").replaceAll("_T0.", ""));
		}
		
		conn.selectOrganism(maize);
		for (Frame f : conn.getAllGFPInstances(geneCategory)) {
			Gene gene = (Gene)f;
//			maizeGeneIDs.add(gene.getLocalID());
//			maizeGeneIDs.add(gene.getCommonName());
			maizeGeneIDs.add(gene.getCommonName().replaceAll("_P0.", "").replaceAll("_T0.", ""));
		}
		
		ArrayList<String> cornMatch = new ArrayList<String>();
		ArrayList<String> cornNoMatch = new ArrayList<String>();
		ArrayList<String> maizeMatch = new ArrayList<String>();
		ArrayList<String> maizeNoMatch = new ArrayList<String>();
		
		conn.selectOrganism(maize);
		for (String id : cornGeneIDs) {
			ArrayList<Frame> searchHits = conn.search(id, geneCategory);
			if (searchHits.isEmpty()) cornNoMatch.add(id);
			else cornMatch.add(id);
		}
		
		conn.selectOrganism(corn);
		for (String id : maizeGeneIDs) {
			ArrayList<Frame> searchHits = conn.search(id, geneCategory);
			if (searchHits.isEmpty()) maizeNoMatch.add(id);
			else maizeMatch.add(id);
		}
		
		System.out.println("Corn genes that match: " + cornMatch.size());
		for (String name : cornMatch) {
			System.out.println(name);
		}
		System.out.println("\n\n");
		
		System.out.println("Corn genes that don't match: " + cornNoMatch.size());
		for (String name : cornNoMatch) {
			System.out.println(name);
		}
		System.out.println("\n\n");
		
		System.out.println("Maize genes that match: " + maizeMatch.size());
		for (String name : maizeMatch) {
			System.out.println(name);
		}
		System.out.println("\n\n");
		
		System.out.println("Maize genes that don't match: " + maizeNoMatch.size());
		for (String name : maizeNoMatch) {
			System.out.println(name);
		}
		System.out.println("\n\n");
	}
	
	static private void geneNameUniqueness(String organism) throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port, "you", "nopass");
		conn.selectOrganism(organism);
		TreeSet<String> geneCommonName = new TreeSet<String>();
		
		ArrayList<Frame> geneFrames = conn.getAllGFPInstances(geneCategory);
		for (Frame f : geneFrames) {
			Gene gene = (Gene)f;
			geneCommonName.add(gene.getCommonName());
//			geneCommonName.add(gene.getCommonName().replaceAll("_P0.", "").replaceAll("_T0.", ""));
		}
		
		System.out.println("Total gene frames: " + geneFrames.size());
		System.out.println("Total unique gene common names: " + geneCommonName.size());
		int i = geneFrames.size() - geneCommonName.size();
		System.out.println("Duplicated common names: " + i);
	}
	
	static private void proteinNameUniqueness(String organism) throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port, "you", "nopass");
		conn.selectOrganism(organism);
		TreeSet<String> proteinCommonName = new TreeSet<String>();
		
		ArrayList<Frame> proteinFrames = conn.getAllGFPInstances(proteinCategory);
		
		System.out.println("Total protein frames: " + proteinFrames.size());
		
		for (Frame f : proteinFrames) {
			Protein protein = (Protein)f;
			proteinCommonName.add(protein.getCommonName());
//			proteinCommonName.add(protein.getCommonName().replaceAll("_P0.", "").replaceAll("_T0.", ""));
		}
		
//		System.out.println("Total protein frames: " + proteinFrames.size());
		System.out.println("Total unique protein common names: " + proteinCommonName.size());
		int i = proteinFrames.size() - proteinCommonName.size();
		System.out.println("Duplicated common names: " + i);
	}
	
	static private void proteinNames(String organism) throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port, "you", "nopass");
		conn.selectOrganism(organism);
		
		String output = "";
		
		@SuppressWarnings("unchecked")
		ArrayList<String> proteinIDs = conn.getClassAllInstances("|Proteins|");
		for (String proteinID : proteinIDs) {
			output += Frame.load(conn, proteinID).getCommonName() + "\n";
		}
		System.out.println(output);
	}
	
	static public ArrayList<String> readFile(String fileName, int column) {
		ArrayList<String> dataRows = new ArrayList<String>();
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(fileName));
			String text = null;
			while ((text = reader.readLine()) != null) {
				try {
					String data = text.split("\t")[column];
					dataRows.add(data);
				} catch (Exception e) {
					dataRows.add("error, column index out of bounds");
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
		
		return dataRows;
	}
}
