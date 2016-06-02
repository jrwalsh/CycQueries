package edu.iastate.cycqueries.queryScripts.DoCornGenesHaveTsAndPs;

import java.util.ArrayList;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.iastate.javacyco.Frame;
import edu.iastate.javacyco.JavacycConnection;
import edu.iastate.javacyco.Network;
import edu.iastate.javacyco.PtoolsErrorException;

public class DoCornGenesHaveTsAndPs {
	private static String hostCorn = "jrwalsh.student.iastate.edu";
	private static String hostMaize = "jrwalsh-server.student.iastate.edu";
	private static String organismCorn = "CORN";
	private static String organismMaize = "MAIZE";
	private static int port = 4444;
	private static String frameType = "|All-Genes|";

	public static void main(String[] args) {
		try {
			Long start = System.currentTimeMillis();
			
//			testTsAndPs(hostCorn, organismCorn);
			testTsAndPs(hostMaize, organismMaize);
			
			Long stop = System.currentTimeMillis();
			Long runtime = (stop - start) / 1000;
			System.out.println("Runtime is " + runtime + " seconds.");
		}
		catch(Exception e) {
			e.printStackTrace();
			System.out.println("Caught a "+e.getClass().getName()+". Shutting down...");
		}
		
	}
	
	private static void testTsAndPs(String host, String organism) throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port);
		System.out.println("Do any genes have both a _P## and _T##, and if so, are they for the same splice variant? Organism = " + organism);
		conn.selectOrganism(organism);
		
		int onlyP = 0;
		int onlyT = 0;
		int both = 0;
		int neither = 0;
		int overlap = 0;
		
		String p1 = "(.*)(_P)(\\d\\d$)";// remove _P##
		Pattern r1 = Pattern.compile(p1);
		
		String p2 = "(.*)(_T)(\\d\\d$)";// remove _T##
		Pattern r2 = Pattern.compile(p2);
		
		String p3 = "(.*)(_FGP)(\\d\\d\\d$)";// remove _FGP###
		Pattern r3 = Pattern.compile(p3);
		
		String p4 = "(.*)(_FGT)(\\d\\d\\d$)";// remove _FGT###
		Pattern r4 = Pattern.compile(p4);
		
		for (String geneName : getAllGeneNames(conn)) {
			boolean hasP = false;
			boolean hasT = false;
			boolean hasOverlap = false;//if the splice number on both a T and a P match, report it
			ArrayList<String> spliceNumbers = new ArrayList<String>();
			ArrayList<Frame> transcripts = conn.search(geneName, frameType);
			for (Frame transcript : transcripts) {
				Matcher m1 = r1.matcher(transcript.getCommonName());
				Matcher m2 = r2.matcher(transcript.getCommonName());
				Matcher m3 = r3.matcher(transcript.getCommonName());
				Matcher m4 = r4.matcher(transcript.getCommonName());
				if (m1.find()) {
					String spliceNumber = m1.group(3);
					if (spliceNumbers.contains(spliceNumber)) hasOverlap = true;
					else spliceNumbers.add(spliceNumber);
					hasP = true;
				} else if (m3.find()) {
					String spliceNumber = m3.group(3);
					if (spliceNumbers.contains(spliceNumber)) hasOverlap = true;
					else spliceNumbers.add(spliceNumber);
					hasP = true;
				} else if (m2.find()) {
					String spliceNumber = m2.group(3);
					if (spliceNumbers.contains(spliceNumber)) hasOverlap = true;
					else spliceNumbers.add(spliceNumber);
					hasT = true;
				} else if (m4.find()) {
					String spliceNumber = m4.group(3);
					if (spliceNumbers.contains(spliceNumber)) hasOverlap = true;
					else spliceNumbers.add(spliceNumber);
					hasT = true;
				}
			}
			if (hasP && hasT) {
				both++;
//				System.out.println(geneName);
//				for (Frame transcript : transcripts) {
//					System.out.println("\t" + transcript.getCommonName());
//				}
			}
			else if (hasP)
				onlyP++;
			else if (hasT)
				onlyT++;
			else
				neither++;
			
			if (hasOverlap) {
				overlap++;
				System.out.println(geneName + " has overlap");
			}
//			System.out.println("Processing : " + geneName);
//			System.out.println(onlyP + "\t" + onlyT + "\t" + both + "\t" + neither);
		}
		System.out.println("OnlyPs\tOnlyTs\tBoth\tNeither\tHasOverlap");
		System.out.println(onlyP + "\t" + onlyT + "\t" + both + "\t" + neither + "\t" + overlap);
	}

	private static TreeSet<String> getAllGeneNames(JavacycConnection conn) throws PtoolsErrorException {
		Network hierarchy = conn.getClassHierarchy(frameType, true, true);
		Set<Frame> allFrames = hierarchy.getNodes();
		
		TreeSet<String> geneNames = new TreeSet<String>();
		
		for (Frame frame : allFrames) {
			String geneName = frame.getCommonName();
			
			if (frame.isClassFrame()) {
//				System.err.println("Ignoring class frame " + frame.getLocalID());
			} else {
				String p1 = "(.*)(_P)(\\d\\d$)";// remove _P##
				Pattern r1 = Pattern.compile(p1);
				Matcher m1 = r1.matcher(geneName);
				m1.find();
				geneName = m1.replaceFirst("$1");
				
				String p2 = "(.*)(_T)(\\d\\d$)";// remove _T##
				Pattern r2 = Pattern.compile(p2);
				Matcher m2 = r2.matcher(geneName);
				m2.find();
				geneName = m2.replaceFirst("$1");
				
				String p3 = "(.*)(_FGP)(\\d\\d\\d$)";// remove _FGP###
				Pattern r3 = Pattern.compile(p3);
				Matcher m3 = r3.matcher(geneName);
				m3.find();
				geneName = m3.replaceFirst("$1");
				
				String p4 = "(.*)(_FGT)(\\d\\d\\d$)";// remove _FGT###
				Pattern r4 = Pattern.compile(p4);
				Matcher m4 = r4.matcher(geneName);
				m4.find();
				geneName = m4.replaceFirst("$1");
				
				geneNames.add(geneName);
			}
		}
		return geneNames;
	}
	
//	public static void test() {
//		JavacycConnection conn = new JavacycConnection(hostCorn, port);
//		System.out.println("Do any genes have both a _P## and _T##, and if so, are they for the same splice variant? Organism = " + organism);
//		conn.selectOrganism(organismCorn);
//		
//		int onlyP = 0;
//		int onlyT = 0;
//		int both = 0;
//		int neither = 0;
//		
//		Network hierarchy = conn.getClassHierarchy(frameType, true, true);
//		Set<Frame> allFrames = hierarchy.getNodes();
//		
//		for (Frame frame : allFrames) {
//			String geneName = frame.getCommonName();
//			
//			if (frame.isClassFrame()) {
//				System.err.println("Ignoring class frame " + frame.getLocalID());
//			} else {
////				System.out.println("Instance frame " + frame.getCommonName());
//				String p1 = "(.*)(_P)(\\d\\d$)";// remove _P##
//				Pattern r1 = Pattern.compile(p1);
//				Matcher m1 = r1.matcher(geneName);
//				m1.find();
//				geneName = m1.replaceFirst("$1");
//				
//				String p2 = "(.*)(_T)(\\d\\d$)";// remove _T##
//				Pattern r2 = Pattern.compile(p2);
//				Matcher m2 = r2.matcher(geneName);
//				m2.find();
//				geneName = m2.replaceFirst("$1");
//				
//				String p3 = "(.*)(_FGP)(\\d\\d\\d$)";// remove _FGP###
//				Pattern r3 = Pattern.compile(p3);
//				Matcher m3 = r3.matcher(geneName);
//				m3.find();
//				geneName = m3.replaceFirst("$1");
//				
//				String p4 = "(.*)(_FGT)(\\d\\d\\d$)";// remove _FGT###
//				Pattern r4 = Pattern.compile(p4);
//				Matcher m4 = r4.matcher(geneName);
//				m4.find();
//				geneName = m4.replaceFirst("$1");
//				
//				ArrayList<Frame> transcriptList = conn.search(geneName, frameType);
//				boolean hasP = false;
//				boolean hasT = false;
//				for (Frame transcript : transcriptList) {
//					m1 = r1.matcher(transcript.getCommonName());
//					m2 = r2.matcher(transcript.getCommonName());
//					m3 = r3.matcher(transcript.getCommonName());
//					m4 = r4.matcher(transcript.getCommonName());
//					if (m1.find() || m3.find()) {
//						hasP = true;
//					}
//					if (m2.find() || m3.find()) {
//						hasT = true;
//					}
//				}
//				
//				if (hasP && hasT)
//					both++;
//				else if (hasP)
//					onlyP++;
//				else if (hasT)
//					onlyT++;
//				else
//					neither++;
//			}
//			System.out.println("Processing : " + frame.getCommonName() + " converted to " + geneName);
//			System.out.println(onlyP + "\t" + onlyT + "\t" + both + "\t" + neither);
//		}
//	}
}
