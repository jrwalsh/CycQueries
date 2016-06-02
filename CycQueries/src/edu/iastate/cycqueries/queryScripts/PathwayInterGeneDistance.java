package edu.iastate.cycqueries.queryScripts;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;

import edu.iastate.javacyco.Frame;
import edu.iastate.javacyco.Gene;
import edu.iastate.javacyco.JavacycConnection;
import edu.iastate.javacyco.Pathway;
import edu.iastate.javacyco.PtoolsErrorException;

@SuppressWarnings("unchecked")
public class PathwayInterGeneDistance {
	private static String host = "localhost";
	private static String organism = "ECOLI";
	private static int port = 4444;
	
	public static void main(String[] args) {
		try {
			Long start = System.currentTimeMillis();
			
			pathwayInterGeneDistance(organism);
			
			Long stop = System.currentTimeMillis();
			Long runtime = (stop - start) / 1000;
			System.out.println("Runtime is " + runtime + " seconds.");
		}
		catch(Exception e) {
			e.printStackTrace();
			System.out.println("Caught a "+e.getClass().getName()+". Shutting down...");
		}
	}
	
	static private void pathwayInterGeneDistance(String organism) throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port);
		conn.selectOrganism(organism);
		
		ArrayList<String> pathwayLabels = (ArrayList<String>) conn.allPathways();
		for (String pathwayLabel : pathwayLabels) {
			Pathway pathway = (Pathway) Pathway.load(conn, pathwayLabel);
			ArrayList<Frame> geneFrames = pathway.getGenes();
			ArrayList<String> positions = new ArrayList<String>();
			int average = getAverage(genesToSortedPositionArray(geneFrames));
			System.out.println(pathway.getCommonName() + "\t" + geneFrames.size() + "\t" + average);
		}
	}
	
	// Convert a list of genes to a list of inter-gene distance values
	static private int[] genesToSortedPositionArray(ArrayList<Frame> genes) {
		int[] positions = new int[genes.size()];
		int index = 0;
		for (Frame gene : genes) {
			try {
				positions[index] = (Integer.parseInt(gene.getSlotValue("LEFT-END-POSITION")));
			} catch (NumberFormatException e) {
				e.printStackTrace();
			} catch (PtoolsErrorException e) {
				e.printStackTrace();
			}
			index++;
		}
		
		Arrays.sort(positions);
		
		int[] diffArray = new int[positions.length-1];
		for (int i = 0; i < positions.length-1; i++) {
			diffArray[i] = positions[i+1] - positions[i];
		}
		return diffArray;
	}
	
	static private int getAverage(int[] values) {
		if (values.length < 1) return 0;
		
		int total = 0;
		for (int i : values) {
			total += i;
		}
		int average = total / values.length; //TODO integer rounding error?
		return average;
	}
}
