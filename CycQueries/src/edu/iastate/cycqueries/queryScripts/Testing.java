package edu.iastate.cycqueries.queryScripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JFileChooser;
import javax.swing.UIManager;

import edu.iastate.javacyco.Frame;
import edu.iastate.javacyco.Gene;
import edu.iastate.javacyco.JavacycConnection;
import edu.iastate.javacyco.Network;
import edu.iastate.javacyco.Protein;
import edu.iastate.javacyco.PtoolsErrorException;
import edu.iastate.javacyco.Reaction;

public class Testing {
	private static String host = "jrwalsh";
	private static String organism = "CORN";
	private static int port = 4444;
	private static JavacycConnection conn;
	private static File logFile;
	
	public static void main(String[] args) {
		
		mergeList();
		
//		try {
//			conn = new JavacycConnection(host, port);
//			conn.selectOrganism(organism);
//			test();
//		} catch (PtoolsErrorException e) {
//			e.printStackTrace();
//		}
	}
	
	public static void mergeList() {
		BufferedReader reader = null;
		String line = null;
		HashMap<String,TreeSet<String>> reviewed = new HashMap<String,TreeSet<String>>(); // UniProt Entry maps to set of EC numbers
		HashMap<String,TreeSet<String>> unreviewed = new HashMap<String,TreeSet<String>>();
		
		try {
			UIManager.setLookAndFeel("com.sun.java.swing.plaf.windows.WindowsLookAndFeel");
			JFileChooser fc = new JFileChooser();
			
			//Process UniProt
			if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) {
				File uniprotFile = fc.getSelectedFile();
				
				reader = new BufferedReader(new FileReader(uniprotFile));
				line = reader.readLine();//skip header
				while ((line = reader.readLine()) != null) {
					String[] data = line.split("\t");
					String entry = data[0];
					String ecNumbers = data[1];
					String existence = data[2];
					String organism = data[3];
					String crefBrenda = data[4];
					String status = data[5];
//					String pubMed = data[6];
					
					if (status.equalsIgnoreCase("reviewed")) {
						if (!reviewed.containsKey(entry)) reviewed.put(entry, new TreeSet<String>());
						for (String ecNumber : ecNumbers.split("; ")) {
							reviewed.get(entry).add(ecNumber); //Keep SwissProt ECs
						}
					} else {
						if (!unreviewed.containsKey(entry)) unreviewed.put(entry, new TreeSet<String>());
//						for (String ecNumber : ecNumbers.split("; ")) {
//							unreviewed.get(entry).add(ecNumber); //Don't keep TrEMBL ECs
//						}
					}
				}
			} else {
				System.err.println("User Canceled");
				return;
			}
			
			//Process BRENDA
			if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) {
				File brendaFile = fc.getSelectedFile();
				
				reader = new BufferedReader(new FileReader(brendaFile));
				line = reader.readLine();//skip header
				while ((line = reader.readLine()) != null) {
					String[] data = line.split("\t");
					String ecNumber = data[0];
					String name = data[1];
					String entry = data[2];
					String organism = data[3];
					
					if (reviewed.containsKey(entry)) {
						reviewed.get(entry).add(ecNumber);
					} else if (unreviewed.containsKey(entry)) {
						unreviewed.get(entry).add(ecNumber);
					} else {
						System.err.println("Unmatched BRENDA entry");// Unmatched values shouldn't happen
					}
					
				}
			} else {
				System.err.println("User Canceled");
				return;
			}
			
			//Print out merged results
			if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION) {
				File outFile = fc.getSelectedFile();
				String outString = "";
				for (String key : reviewed.keySet()) {
					outString += key + "\t" + reviewed.get(key) + "\treviewed\n";
				}
				for (String key : unreviewed.keySet()) {
					outString += key + "\t" + unreviewed.get(key) + "\tunreviewed\n";
				}
				appendLine(outFile, outString);
			} else {
				System.err.println("User Canceled");
				return;
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
	
	public static void test() throws PtoolsErrorException {
		logFile = new File("C:\\Users\\Jesse\\Desktop\\New folder", "out.tab");
//		ArrayList<Frame> frames = conn.search("GRMZM2G119578_P01", "|All-Genes|");
//		frames.get(0).print();
		
//		Frame.load(conn, "GDQC-109058-MONOMER").print();
		
		Network hierarchy = conn.getClassHierarchy("|Genes|", true, true);
		Set<Frame> allFrames = hierarchy.getNodes();
		for (Frame frame : allFrames) {
			if (!frame.isClassFrame()) {
				for (String EC : getECList(conn, (Gene)frame)) {
					appendLine(logFile, frame.getLocalID() + "\t" + frame.getCommonName() + "\t" + EC + "\n");
				}
			}
		}
		
	}
	
	private static TreeSet<String> getECList(JavacycConnection conn, Gene gene) throws PtoolsErrorException {
		TreeSet<String> ECList = new TreeSet<String>();
		ArrayList<Protein> proteins = gene.getProducts();//Genes can have multiple proteins, and rarely no proteins or non-protein products
		if (proteins.size() == 0) {
			System.err.println("No protein");
			return null;
		} else {
			for (Protein protein : proteins) {
				for (Frame catalysis : protein.getCatalysis()) {
					boolean computationalEVCode = false;
					ArrayList<String> citations = catalysis.getSlotValues("CITATIONS");
					if (citations == null || citations.isEmpty()) computationalEVCode = true; //MaizeCyc is all computational, and does not included it explicitly
					else {
						for (String citation : citations) {
							if (citation.toString().contains("EV-COMP")) { //We don't care how many citations, as long as one of them indicates the EC was computationally predicted
								computationalEVCode = true;
							}
						}
					}
					if (computationalEVCode) {
						Reaction reaction = (Reaction) Reaction.load(conn, catalysis.getSlotValue("REACTION"));
						String EC = null;
						if (reaction != null) EC = reaction.getEC();
						if (EC != null && EC.length() > 0) {
							EC = EC.replace("EC-", "");//Remove the EC- prefix
							String p = "\\d{1,3}?.\\d{1,3}?.\\d{1,3}?.\\d{1,3}?";//EC numbers are, like IP addresses, 4 numbers of up to 3 digits each, separated by periods. Don't allow partial EC numbers, such as 1.-.-.- or 1.2.3.-
							Pattern r = Pattern.compile(p);
							Matcher m = r.matcher(EC);
							if (m.find()) {
								ECList.add(EC);//Keep full EC annotations if there was either some computational evidence or no evidence at all (assume computational)
							} else {
								System.err.println("partial EC? : " + EC);
							}
							
						}
					} else {
						System.err.println("Failed to find computational citation for :" + gene.getLocalID());
					}
				}
			}
		}
		return ECList;
	}
	
//	private static TreeSet<String> getECList(JavacycConnection conn, Gene gene) throws PtoolsErrorException {
//		TreeSet<String> ECList = new TreeSet<String>();
//		ArrayList<Protein> proteins = gene.getProducts();//Genes can have multiple proteins, and rarely no proteins or non-protein products
//		if (proteins.size() == 0) {
//			System.err.println("No protein");
//			return null;
//		} else {
//			for (Protein protein : proteins) {
//				for (Frame catalysis : protein.getCatalysis()) {
//					boolean computationalEVCode = false;
//					ArrayList<String> citations = catalysis.getSlotValues("CITATIONS");
//					if (citations == null || citations.isEmpty()) computationalEVCode = true; //MaizeCyc is all computational, and does not included it explicitly
//					else {
//						for (String citation : citations) {
//							if (citation.toString().contains("EV-COMP")) { //We don't care how many citations, as long as one of them indicates the EC was computationally predicted
//								computationalEVCode = true;
//							}
//						}
//					}
//					if (computationalEVCode) {
//						Reaction reaction = (Reaction) Reaction.load(conn, catalysis.getSlotValue("REACTION"));
//						String EC = null;
//						if (reaction != null) EC = reaction.getEC();
//						if (EC != null && EC.length() > 0) {
//							EC = EC.replace("EC-", "");//Remove the EC- prefix
//							String p = "\\d{1,3}?.\\d{1,3}?.\\d{1,3}?.\\d{1,3}?";//EC numbers are, like IP addresses, 4 numbers of up to 3 digits each, separated by periods. Don't allow partial EC numbers, such as 1.-.-.- or 1.2.3.-
//							Pattern r = Pattern.compile(p);
//							Matcher m = r.matcher(EC);
//							if (m.find()) {
//								ECList.add(EC);//Keep full EC annotations if there was either some computational evidence or no evidence at all (assume computational)
//							} else {
//								System.err.println("partial EC? : " + EC);
//							}
//							
//						}
//					} else {
//						System.err.println("Failed to find computational citation for :" + gene.getLocalID());
//					}
//				}
//			}
//		}
//		return ECList;
//	}
	
	protected static void appendLine(File file, String printString) {
		PrintStream o = null;
		try {
			o = new PrintStream(new FileOutputStream(file, true));
			o.append(printString);
			o.close();
		}
		catch(Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
	}
}
