package edu.iastate.cycqueries.queryScripts;

import java.util.ArrayList;

import edu.iastate.javacyco.Frame;
import edu.iastate.javacyco.JavacycConnection;
import edu.iastate.javacyco.PtoolsErrorException;

public class Testing {
	private static String host = "jrwalsh-server.student.iastate.edu";
	private static String organism = "MAIZE";
	private static int port = 4444;
	private static JavacycConnection conn;
	
	public static void main(String[] args) {
		conn = new JavacycConnection(host, port);
		conn.selectOrganism(organism);
		
		try {
			test();
		} catch (PtoolsErrorException e) {
			e.printStackTrace();
		}
	}
	
	public static void test() throws PtoolsErrorException {
//		ArrayList<Frame> frames = conn.search("GRMZM2G119578_P01", "|All-Genes|");
//		frames.get(0).print();
		
		Frame.load(conn, "ENZRXNBWI-261").print();
	}
}
