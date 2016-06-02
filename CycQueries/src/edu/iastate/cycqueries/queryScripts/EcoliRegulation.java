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
public class EcoliRegulation {
	private static String host = "jrwalsh.student.iastate.edu";
	private static String organism = "ECOLI";
	private static int port = 4444;
	
	public static void main(String[] args) {
		try {
			Long start = System.currentTimeMillis();
			
			regulationQuery();
			
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
	static private void regulationQuery() throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port);
		conn.selectOrganism(organism);
		
		ArrayList<Frame> regulators = new ArrayList<Frame>();
		regulators.add(Frame.load(conn, "PD00266"));
		regulators.add(Frame.load(conn, "ALAS-CPLX"));
		regulators.add(Frame.load(conn, "CPLX0-2021"));
		regulators.add(Frame.load(conn, "G7678-MONOMER"));
//		regulators.add(Frame.load(conn, "EG10821"));
//		regulators.add(Frame.load(conn, "EG10935"));
//		regulators.add(Frame.load(conn, "G6420"));
//		regulators.add(Frame.load(conn, "EG10359"));
//		regulators.add(Frame.load(conn, "EG10440"));
		
		for (Frame regulator : regulators) {
			ArrayList<String> regulates = regulator.getSlotValues("Regulates");
			for (String regulate : regulates) {
				Frame.load(conn, regulate).print();
			}
		}
		
	}
	
}
