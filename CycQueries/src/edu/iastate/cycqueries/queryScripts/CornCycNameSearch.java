package edu.iastate.cycqueries.queryScripts;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeSet;
import java.util.Map.Entry;
import java.util.regex.Pattern;

import edu.iastate.javacyco.Frame;
import edu.iastate.javacyco.Gene;
import edu.iastate.javacyco.JavacycConnection;
import edu.iastate.javacyco.Protein;
import edu.iastate.javacyco.PtoolsErrorException;

/**
 * Written for Taner Sen of the MaizeGDB/CornCyc group.
 * 
 * 
 */
public class CornCycNameSearch {
	private static String host = "jrwalsh.student.iastate.edu";
	private static String organismCorn = "CORN";
	private static String geneCategory = "|All-Genes|";
	private static String proteinCategory = "|Proteins|";
	private static int port = 4444;
	private static String fileName = "C:\\Users\\Jesse\\Desktop\\test.txt";
	
	public static void main(String[] args) {
		try {
			Long start = System.currentTimeMillis();
			
			getMatchingGeneFrames_version2();
			
			Long stop = System.currentTimeMillis();
			Long runtime = (stop - start) / 1000;
			System.out.println("Runtime is " + runtime + " seconds.");
		}
		catch(Exception e) {
			e.printStackTrace();
			System.out.println("Caught a "+e.getClass().getName()+". Shutting down...");
		}
	}
	
	static private void getMatchingGeneFrames() throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port);
		conn.selectOrganism(organismCorn);
		
		ArrayList<String> noMatches = new ArrayList<String>();
		HashMap<String, String> singleMatches = new HashMap<String, String>();
		HashMap<String, ArrayList<String>> multipleMatches = new HashMap<String, ArrayList<String>>();
		
		for (String name : readFile(fileName, 0)) {
			ArrayList<Frame> searchHits = conn.search(name, proteinCategory);
			if (searchHits.isEmpty()) {
				noMatches.add(name);
			} else if (searchHits.size() == 1) singleMatches.put(name, searchHits.get(0).getLocalID());
			else {
				ArrayList<String> matches = new ArrayList<String>();
				for (Frame hit : searchHits) matches.add(hit.getLocalID());
				multipleMatches.put(name, matches);
			}
		}
		
		System.out.println("Names with no matchs: " + noMatches.size());
		for (String name : noMatches) System.out.println(name);
		System.out.println("\n\n");
				
		System.out.println("Names with single matchs: " + singleMatches.keySet().size());
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
	
	// Run on 3/21/2014 to fill Taner's request for partial matches in CornCyc to several alternative identifiers used in GO annotations
	static private void getMatchingGeneFrames_version2() throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port);
		conn.selectOrganism(organismCorn);
		
		TreeSet<String> inputIDs = new TreeSet<String>();
		for (String userProvidedID : readFile(fileName, 0)) {
			inputIDs.add(userProvidedID);
		}
		
		// For each user provided ID, do a search of that ID in the current DB.  Store results.
		HashMap<String, ArrayList<Frame>> searchResults = new HashMap<String, ArrayList<Frame>>();
		for (String userProvidedID : inputIDs) {
			try {
				if (!searchResults.containsKey(userProvidedID)) {
					searchResults.put(userProvidedID, conn.search(userProvidedID, proteinCategory));
				}
			} catch (PtoolsErrorException e) {
				e.printStackTrace();
				searchResults.put(userProvidedID, null);
			}
		}
		
		// Sort out the original list into different groups based on the quality of the match
		for (String userProvidedID : inputIDs) {
			if (searchResults.get(userProvidedID) == null || searchResults.get(userProvidedID).size() == 0) {
				System.out.println(userProvidedID + "\tNo Match\t"); // Case: No matches returned by search
			} else if (conn.frameExists(userProvidedID)) {
				System.out.println(userProvidedID + "\tFrameID Match\t" + userProvidedID); // Case: No matches returned by search
			} else {
				int synonymMatchesCount = 0;
				ArrayList<String> exactFrameMatches = new ArrayList<String>();
				ArrayList<String> substringFrameMatches = new ArrayList<String>();
				for (Frame match : searchResults.get(userProvidedID)) {
					try {
						for (String name : match.getNames()) {
							if (userProvidedID.equalsIgnoreCase(name.replaceAll("\"", ""))) {
								synonymMatchesCount++;
								exactFrameMatches.add(match.getLocalID() + "(" + name + ")");
								break; // Only increment by 1 per frame with a matching synonym
							} else if (Pattern.compile(Pattern.quote(userProvidedID), Pattern.CASE_INSENSITIVE).matcher(name).find()) {//(name.contains(userProvidedID)) {
								substringFrameMatches.add(match.getLocalID() + "(" + name + ")");
							}
						}
					} catch (PtoolsErrorException e) {
						e.printStackTrace();
					}
				}
				
				if (synonymMatchesCount == 0) {
					if (substringFrameMatches.size() == 0) {
//						System.err.println("error 1 " + userProvidedID);
						for (Frame match : searchResults.get(userProvidedID)) {
							System.out.println(userProvidedID + "\tSubstring Synonym Match\t" + match.getLocalID()); // Case: No exact matches found by search, but partial matches found. Can't identify which name the partial match matches to
						}
					}
					for (String hit : substringFrameMatches) System.out.println(userProvidedID + "\tSubstring Synonym Match\t" + hit); // Case: No exact matches found by search, but partial matches found
				} else if (synonymMatchesCount == 1) {
//					if (exactFrameMatches.size() != 1) System.err.println("error 2 " + userProvidedID);
					for (String hit : exactFrameMatches) System.out.println(userProvidedID + "\tSingle Exact Synonym Match\t" + hit); // Case: 1 and only 1 synonym had an exact match
				} else {
//					if (exactFrameMatches.size() == 0) System.err.println("error 3 " + userProvidedID);
					for (String hit : exactFrameMatches) System.out.println(userProvidedID + "\tMultiple Exact Synonym Match\t" + hit); // Case: We have exact matches on the synonyms of two or more frames.  Unlikely, but if it happens we can't use this match.
				}
			}
		}
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
