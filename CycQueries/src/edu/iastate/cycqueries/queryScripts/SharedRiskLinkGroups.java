package edu.iastate.cycqueries.queryScripts;

import java.util.ArrayList;
import java.util.HashMap;

import edu.iastate.javacyco.EnzymeReaction;
import edu.iastate.javacyco.Frame;
import edu.iastate.javacyco.Gene;
import edu.iastate.javacyco.JavacycConnection;
import edu.iastate.javacyco.Pathway;
import edu.iastate.javacyco.Protein;
import edu.iastate.javacyco.PtoolsErrorException;
import edu.iastate.javacyco.Reaction;
import edu.iastate.javacyco.TranscriptionUnit;

public class SharedRiskLinkGroups {
	private static String host = "localhost";
	private static String organism = "ECOLI";
	private static String geneCategory = "|All-Genes|";
	private static String transcriptionUnitCategory = "|Transcription-Units|";
	private static String reactionCategory = "|Small-Molecule-Reactions|";//"|Chemical-Reactions|"
	private static int port = 4444;
	
	
	public static void main(String[] args) {
		try {
			Long start = System.currentTimeMillis();
//			genesInTranscriptionUnit();
			averageGeneSpacingInPathway();
			Long stop = System.currentTimeMillis();
			Long runtime = (stop - start) / 1000;
			System.out.println("Runtime is " + runtime + " seconds.");
		}
		catch(Exception e) {
			e.printStackTrace();
			System.out.println("Caught a "+e.getClass().getName()+". Shutting down...");
		}
	}
	
	static private void genesInTranscriptionUnit() throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port, "you", "nopass");
		conn.selectOrganism(organism);
		
		int srlgCount = 0;
		
		ArrayList<Frame> transcriptionUnits = conn.getAllGFPInstances(transcriptionUnitCategory);
		for (Frame transcriptionUnit : transcriptionUnits) {
			TranscriptionUnit tu = (TranscriptionUnit) transcriptionUnit;
			HashMap<String, ArrayList<String>> genesToReactions = new HashMap<String, ArrayList<String>>();
			ArrayList<Gene> genes = tu.getGenes();
			for (Gene gene : genes) {
				ArrayList<Protein> enzymes = gene.getEnzymes();
				ArrayList<String> reactions = new ArrayList<String>();
				for (Protein enzyme : enzymes) {
					reactions.addAll(conn.reactionsOfEnzyme(enzyme.getLocalID()));
				}
				genesToReactions.put(gene.getLocalID(), reactions);
			}
			
			int reactionCount = 0;
			boolean srlgFound = false;
			String printString = "";
			printString += tu.getLocalID() + "\n";
			for (String geneID : genesToReactions.keySet()) {
				printString += "\t" + geneID + "\n";
				
				ArrayList<String> reactionIDs = genesToReactions.get(geneID);
				reactionCount += reactionIDs.size();
				if (reactionCount > 1 && !srlgFound) {
					// Only count once per tu
					srlgFound = true;
					srlgCount++;
				}
				
				for (String reactionID : genesToReactions.get(geneID)) {
					printString += "\t\t" + reactionID + "\n";
				}
			}
			if (srlgFound) System.out.print(printString);
		}
		System.out.println("Found " + srlgCount + " SRLGs.");
	}
	
	static private void averageGeneSpacingInPathway() throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port, "you", "nopass");
		conn.selectOrganism(organism);
		
		ArrayList<String> pathways = conn.allPathways();
		
		for (String pathway : pathways) {
			Pathway pwy = (Pathway) Pathway.load(conn, pathway);
			ArrayList<Frame> geneFrames = pwy.getGenes();
			
			System.out.println(pathway);
			for (Frame geneFrame : geneFrames) {
				Gene g = (Gene)geneFrame;
				System.out.println("\t" + g.getCommonName() + "\t" + g.getSlotValue("LEFT-END-POSITION") + "\t" + g.getSlotValue("RIGHT-END-POSITION"));
			}
		}
	}
	
	static private int distanceBetweenGenes(Gene gene1, Gene gene2) throws NumberFormatException, PtoolsErrorException {
		int gene1left = Integer.parseInt(gene1.getSlotValue("LEFT-END-POSITION"));
		int gene1right = Integer.parseInt(gene1.getSlotValue("RIGHT-END-POSITION"));
		int gene2left = Integer.parseInt(gene2.getSlotValue("LEFT-END-POSITION"));
		int gene2right = Integer.parseInt(gene2.getSlotValue("RIGHT-END-POSITION"));
		
		// For this application, we don't particularly care which direction the gene is transcribed.  We merely want to find the size of space between them.
		// Microbial genomes are circular, so if left position is greater than right, than this gene crosses the beginning of the genome
		// If one gene is on the end, and the other at the beginning, the distance will be huge, when in reality it should be small.  I'll risk it for testing atm.
		
		// If genes overlap, than return 0 distance
		if (gene1left <= gene2right && gene2left <= gene1right)
			return 0;
		// Don't know which gene is first, but don't return a negative distance
		else if (gene1right < gene2left)
			return gene2left-gene1right;
		else 
			//if (gene2right < gene1left)
			return gene1left-gene1right;
	}
	
	static private boolean isAnyGeneContainedByAnother(ArrayList<Gene> genes) throws NumberFormatException, PtoolsErrorException {
		while (genes.size() > 0) {
			Gene thisGene = genes.remove(0);
			int thisGeneStart = Integer.parseInt(thisGene.getSlotValue("LEFT-END-POSITION"));
			int thisGeneEnd = Integer.parseInt(thisGene.getSlotValue("RIGHT-END-POSITION"));
			for (Gene gene : genes) {
				int geneStart = Integer.parseInt(gene.getSlotValue("LEFT-END-POSITION"));
				int geneEnd = Integer.parseInt(gene.getSlotValue("RIGHT-END-POSITION"));
				
				if (containsRange(thisGeneStart, thisGeneEnd, geneStart, geneEnd))
					return true;
			}
		}
		
		return false;
	}
	
	static private boolean containsRange(int start1, int end1, int start2, int end2) {
		if (start1 <= start2 && end1 >= end2)
			return true;
		else if (start2 <= start1 && end2 >= end1)
			return true;
		else
			return false;
	}
	
 	static private void reactionsOfGenes() throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port, "you", "nopass");
		conn.selectOrganism(organism);
		
		ArrayList<Frame> genes = conn.getAllGFPInstances(geneCategory);
		for (Frame gene : genes) {
			Gene g = (Gene) gene;
			ArrayList<Protein> enzymes = g.getEnzymes();
			for (Protein enzyme : enzymes) {
				ArrayList<String> reactions = conn.reactionsOfEnzyme(enzyme.getLocalID());
				for (String reaction : reactions) {
					System.out.println(g.getLocalID() + "\t" + reaction);
				}
			}
			
		}
	}
	
	static private String genesOfReaction(String reactionID, JavacycConnection conn) throws PtoolsErrorException {
		String orRule = "";
		for (Object enzyme : conn.enzymesOfReaction(reactionID)) {
			String andRule = "";
			for (Object gene : conn.genesOfProtein(enzyme.toString())) {
				String geneID = gene.toString();
				try {
					geneID = Frame.load(conn, geneID).getSlotValue("ACCESSION-1").replace("\"", "");
				} catch (Exception e) {
					geneID = gene.toString();
				}
				andRule += geneID + " and ";
			}
			if (andRule.length() > 0) {
				andRule = "(" + andRule.substring(0, andRule.length()-5) + ")";
				orRule += andRule + " or ";
			}
		}
		if (orRule.length() > 0) orRule = orRule.substring(0, orRule.length()-4);
		return orRule;
	}
}
