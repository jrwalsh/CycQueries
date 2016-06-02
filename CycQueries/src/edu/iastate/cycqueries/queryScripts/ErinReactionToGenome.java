package edu.iastate.cycqueries.queryScripts;

import java.util.ArrayList;

import edu.iastate.javacyco.Catalysis;
import edu.iastate.javacyco.Frame;
import edu.iastate.javacyco.Gene;
import edu.iastate.javacyco.JavacycConnection;
import edu.iastate.javacyco.Promoter;
import edu.iastate.javacyco.Protein;
import edu.iastate.javacyco.PtoolsErrorException;
import edu.iastate.javacyco.Reaction;
import edu.iastate.javacyco.Regulation;
import edu.iastate.javacyco.TranscriptionUnit;

public class ErinReactionToGenome {
	private static String host = "jrwalsh.student.iastate.edu";
	private static String organism = "ECOLI";
	private static int port = 4444;
	
	private static String user = "me";
	private static String password = "pass";
	
	public static void main(String[] args) {
		try {
			if(args.length<1) {
				System.out.println("Usage: Main REACTION_ID");
				System.exit(0);
			}
			String reactionID = args[0];
			
//			String reactionID = "PGLUCISOM-RXN";
			
			Long start = System.currentTimeMillis();
			queryReactionToGenome(reactionID);
			Long stop = System.currentTimeMillis();
			Long runtime = (stop - start) / 1000;
			System.out.println("Runtime is " + runtime + " seconds.");
		}
		catch(Exception e) {
			e.printStackTrace();
			System.out.println("Caught a "+e.getClass().getName()+". Shutting down...");
		}
	}
	
	static private void queryReactionToGenome(String reactionID) throws PtoolsErrorException {
//		JavacycConnection conn = new JavacycConnection(host, port);
		JavacycConnection conn = new JavacycConnection(host, port, user, password);
		conn.selectOrganism(organism);
		
		Reaction reaction = (Reaction) Reaction.load(conn, reactionID);
		Catalysis eReaction = (Catalysis) Catalysis.load(conn, reaction.getSlotValue("ENZYMATIC-REACTION"));
		Protein enzyme = eReaction.getEnzyme();
		
		ArrayList<Gene> genes = enzyme.getGenes();
		for (Gene gene : genes) {
			ArrayList<TranscriptionUnit> tus = gene.getTranscriptionUnits();
			
			for (TranscriptionUnit tu : tus) {
				tu.print();
				Promoter promoter = tu.getPromoter();
				ArrayList<Frame> regulations = promoter.getRegulations();
				for (Frame regulation : regulations) {
					Frame regulator = ((Regulation)regulation).getRegulator();
					regulator.print();
				}
			}
		}
	}
}
