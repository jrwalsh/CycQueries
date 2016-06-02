package edu.iastate.cycqueries.queryScripts;

import java.util.ArrayList;
import edu.iastate.javacyco.Frame;
import edu.iastate.javacyco.JavacycConnection;
import edu.iastate.javacyco.Pathway;
import edu.iastate.javacyco.PtoolsErrorException;


/**
 * Written for Ruolin.
 * 
 * Designed to query all genes in pathways Ruolin associates with chlorophyll.
 * Returns a tab delimited file with PathwayFrameID, GeneFrameID
 */
public class RuolinChlorophyllGenes {
	private static String host = "localhost";
	private static String organism = "MAIZE";
//	private static String organism = "CORN";
	private static int port = 4444;
	private static String frameType = "|Chlorophyll-Biosynthesis|";
	
	public static void main(String[] args) {
		try {
			Long start = System.currentTimeMillis();
			queryChlorophyllGenes();
			Long stop = System.currentTimeMillis();
			Long runtime = (stop - start) / 1000;
			System.out.println("Runtime is " + runtime + " seconds.");
		}
		catch(Exception e) {
			e.printStackTrace();
			System.out.println("Caught a "+e.getClass().getName()+". Shutting down...");
		}
	}
	
	static private void queryChlorophyllGenes() throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port);
		conn.selectOrganism(organism);
		
		ArrayList<Frame> pahways = getAllOfFrameType(frameType, conn);
		pahways.add(Pathway.load(conn, "PWY-5098"));
		
		for (Frame pathway : pahways) {
			String genes = genesOfPathwayTabDelimited((Pathway) pathway, conn);
			
			System.out.println(pathway.getLocalID() + "\t" + pathway.getCommonName() + "\t" + genes);
		}
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
	
	static private String genesOfPathwayTabDelimited(Pathway pathway, JavacycConnection conn) {
		String returnString = "";
		try {
			ArrayList<Frame> genes = pathway.getGenes();
			for (int i = 0; i < genes.size(); i++) {
				returnString += genes.get(i).getCommonName();
				if (i < genes.size()) returnString += "\t";
			}
		} catch (PtoolsErrorException e) {
			e.printStackTrace();
		}
		return returnString;
	}
}
