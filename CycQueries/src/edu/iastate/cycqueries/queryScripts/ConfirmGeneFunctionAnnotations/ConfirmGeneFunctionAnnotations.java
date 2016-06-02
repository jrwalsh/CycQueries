package edu.iastate.cycqueries.queryScripts.ConfirmGeneFunctionAnnotations;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.iastate.javacyco.Frame;
import edu.iastate.javacyco.Gene;
import edu.iastate.javacyco.JavacycConnection;
import edu.iastate.javacyco.Protein;
import edu.iastate.javacyco.PtoolsErrorException;
import edu.iastate.javacyco.Reaction;

/**
 * Given a set of maize genes and their putative EC values obtained by blasting manually reviewed, experimentally determined maize proteins in UniProt against the B73_v2 maize
 * reference sequence, we would like to confirm that these EC values were computationally assigned in CornCyc and/or MaizeCyc such that we can calculate the TP/FP/TN/FN values
 * for each database as described in the corncompare paper.
 * This is to be done on the publically available base versions of each database, not including any schema upgrades or update propagation steps.
 * -- (i.e. CornCyc 4.0 in 16.5 and MaizeCyc 2.2 in 17.0)
 * 
 * !! The search criteria for MaizeCyc differs from CornCyc, as MaizeCyc does not include accurate splice information, the _P## will not match, and also
 * MaizeCyc does not have any citation data in the Enzymatic Reactions.  The code was manually modified to account for these two discrepancies.
 * 
 * @throws PtoolsErrorException
 * @author Jesse Walsh 5/31/2016
 */
public class ConfirmGeneFunctionAnnotations {
	private static String hostCorn = "jrwalsh.student.iastate.edu";
	private static String hostMaize = "jrwalsh-server.student.iastate.edu";
	private static String organismCorn = "CORN";
	private static String organismMaize = "MAIZE";
	private static int port = 4444;
	private static String frameType = "|All-Genes|";
	private static String inFile = "C:\\Users\\Jesse\\workspace\\CycQueries\\bin\\edu\\iastate\\cycqueries\\queryScripts\\ConfirmGeneFunctionAnnotations\\blast-results-fixed";
	//We expect this inFile to be in the format 
	//Protein \t EC-Numbers
	//Where protein is the blast hit result which represents a corn gene model name and is used as the query term
	//and EC-Numbers is a semi-colon delimited list of EC numbers uniprot associates with this protein
	
	public static void main(String[] args) {
		try {
			Long start = System.currentTimeMillis();
			
			listEVCodes(hostCorn, organismCorn);
//			listEVCodes(hostMaize, organismMaize);
			
//			listEVCodesAllowMultipleGeneHits(hostCorn, organismCorn);
			
			Long stop = System.currentTimeMillis();
			Long runtime = (stop - start) / 1000;
			System.out.println("Runtime is " + runtime + " seconds.");
		}
		catch(Exception e) {
			e.printStackTrace();
			System.out.println("Caught a "+e.getClass().getName()+". Shutting down...");
		}
		
	}
	
	/**
	 * Take in a list of gene model name to EC pairs (with additional columns of info) and
	 * 1) search gene model name to resolve into a gene
	 * 2) get products of gene
	 * 3) get catalyzes of product (i.e. enzymatic reactions)
	 * 4) print citations for the enzymatic reactions
	 * 5) get reaction of enzymatic reaction
	 * 6) print EC number of reaction
	 * 
	 * @throws PtoolsErrorException
	 */
	private static void listEVCodes(String host, String organism) throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port);
		System.out.println("Confirm EV codes for Blast Results: " + organism);
		System.out.println("ProteinFrameID\tEC-Number");
		conn.selectOrganism(organism);
		
		int TP = 0;//When the uniprot ec matches the cyc database EC
		int FP = 0;//When the cyc database has an EC that is not in the Uniprot EC list
		int FN = 0;//When Uniprot has an EC that is not in the cyc database EC list
		int TN = 0;//We cannot capture this value with the data we have
		
		BufferedReader reader = null;
		String line = null;
		try {
			reader = new BufferedReader(new FileReader(inFile));
			reader.readLine();//Clear header line
			
			while ((line = reader.readLine()) != null) {
				String[] data = line.split("\t");
				String item = data[0];
				TreeSet<String> UniProtECs = new TreeSet<String>();
				for (String EC : data[1].split("; ")) {
					UniProtECs.add(EC);
				}

				Gene gene = getGeneMatch(conn, item);//Find the gene which matches the uniprot blast hit
				if (gene == null) {
					System.err.println("No gene match for " + item);
					continue;
				} else {
					TreeSet<String> ECList = getECList(conn, gene);

					//Print the list of EC values for the cyc database
					if (ECList == null || ECList.isEmpty()) {
						System.out.println(gene.getLocalID() + "\t");
					} else {
						for (String EC : ECList) {
							System.out.println(gene.getLocalID() + "\t" + EC);
						}
					}
					
					//Calculate TP,FP,FN
					for (String uniprotEC : UniProtECs) {
						if (ECList.contains(uniprotEC)) TP++;
						else FN++;
					}
					for (String predictedEC : ECList) {
						if (UniProtECs.contains(predictedEC)) {
							//ignore, already counted TPs
						} else FP++;
					}
				}
			}
			
			System.out.println("TP = " + TP);
			System.out.println("FP = " + FP);
			System.out.println("FN = " + FN);
			System.out.println("TN cannot be determined");
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
	
	private static Gene getGeneMatch(JavacycConnection conn, String query) throws PtoolsErrorException {
		ArrayList<Frame> list = conn.search(query, frameType); // Search for gene
		if (list.size() == 1) {
			return (Gene) list.get(0);
		} else if (list.size() == 0) {
			//Try to replace any _P## with an _T##, common error when searching CornCyc/MaizeCyc
			String p = "(.*)(_P)(\\d\\d$)"; // 3 capture groups.  first group is any length of letters or numbers, followed by an _P, followed by 2 numbers and an end of line
			Pattern r = Pattern.compile(p);
			Matcher m = r.matcher(query);
			if (m.find()) {
				query = m.replaceFirst("$1_T$3");
//				query = m.replaceFirst("$1");//MaizeCyc does not differentiate between _P## variants
			} else {
				//The other common reason for not finding a search term: Replace any _FGP### with an _FGT###
				p = "(.*)(_FGP)(\\d\\d\\d$)";
				r = Pattern.compile(p);
				m = r.matcher(query);
				if (m.find()) {
					query = m.replaceFirst("$1_FGT$3");
				}
			}
			list = conn.search(query, frameType); // Search for gene using modified query
			if (list.size() == 1) {
				return (Gene) list.get(0);
			}
		}
		return null; //Search fails if it finds more than 1 hit on initial search, or if it finds 0 or more than 1 hit on modified search
	}
	
	/**
	 * Similar to above method, but if the search returns multiple hits (i.e. splice variants), merge them all together (i.e. do comparison at the gene level, not the
	 * splice variant level).  To do this, remove the _P##, _T##, _FGP####, and _FGT####.
	 * @param host
	 * @param organism
	 * @throws PtoolsErrorException
	 */
	private static void listEVCodesAllowMultipleGeneHits(String host, String organism) throws PtoolsErrorException {
		JavacycConnection conn = new JavacycConnection(host, port);
		System.out.println("Confirm EV codes for Blast Results, merge splice variants: " + organism);
		System.out.println("ProteinFrameID\tEC-Number");
		conn.selectOrganism(organism);
		
		int TP = 0;//When the uniprot ec matches the cyc database EC
		int FP = 0;//When the cyc database has an EC that is not in the Uniprot EC list
		int FN = 0;//When Uniprot has an EC that is not in the cyc database EC list
		int TN = 0;//We cannot capture this value with the data we have
		
		BufferedReader reader = null;
		String line = null;
		try {
			reader = new BufferedReader(new FileReader(inFile));
			reader.readLine();//Clear header line
			
			while ((line = reader.readLine()) != null) {
				String[] data = line.split("\t");
				String item = data[0];
				TreeSet<String> UniProtECs = new TreeSet<String>();
				for (String EC : data[1].split("; ")) {
					UniProtECs.add(EC);
				}

				ArrayList<Gene> genes = getGeneMatchMultiple(conn, item);//Find the gene which matches the uniprot blast hit
				if (genes == null || genes.isEmpty()) {
					System.err.println("No gene match for " + item);
					continue;
				} else {
//					System.out.println(item + " matched up with " + genes.size() + " hits");
					TreeSet<String> ECList = new TreeSet<String>();
					for (Gene gene : genes) {
						TreeSet<String> temp = getECList(conn, gene);
						if (temp != null && !temp.isEmpty()) {
							ECList.addAll(temp);
						}
					}

					//Print the list of EC values for the cyc database
					if (ECList == null || ECList.isEmpty()) {
						System.out.println(item + "\t");
					} else {
						for (String EC : ECList) {
							System.out.println(item + "\t" + EC);
						}
					}
					
					//Calculate TP,FP,FN
					for (String uniprotEC : UniProtECs) {
						if (ECList.contains(uniprotEC)) TP++;
						else FN++;
					}
					for (String predictedEC : ECList) {
						if (UniProtECs.contains(predictedEC)) {
							//ignore, already counted TPs
						} else FP++;
					}
				}
			}
			
			System.out.println("TP = " + TP);
			System.out.println("FP = " + FP);
			System.out.println("FN = " + FN);
			System.out.println("TN cannot be determined");
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
	
	private static ArrayList<Gene> getGeneMatchMultiple(JavacycConnection conn, String query) throws PtoolsErrorException {
		ArrayList<Gene> genes = new ArrayList<Gene>();
		
		String p = "(.*)(_P)(\\d\\d$)"; // 3 capture groups.  first group is any length of letters or numbers, followed by an _P, followed by 2 numbers and an end of line
		Pattern r = Pattern.compile(p);
		Matcher m = r.matcher(query);
		if (m.find()) {
			query = m.replaceFirst("$1");
		} else {
			//The other common reason for not finding a search term: Replace any _FGP### with an _FGT###
			p = "(.*)(_FGP)(\\d\\d\\d$)";
			r = Pattern.compile(p);
			m = r.matcher(query);
			if (m.find()) {
				query = m.replaceFirst("$1");
			}
		}
		ArrayList<Frame> list = conn.search(query, frameType); // Search for gene using modified query
		
		if (list.size() > 0) {
			for (Frame frame : list) genes.add((Gene)frame);
			return genes;
		}
		return null;
	}
	
	private static TreeSet<String> getECList(JavacycConnection conn, Gene gene) throws PtoolsErrorException {
		TreeSet<String> ECList = new TreeSet<String>();
		ArrayList<Protein> proteins = gene.getProducts();//Genes can have multiple proteins, and rarely no proteins or non-protein products
		if (proteins.size() == 0) {
			System.err.println("No protein");
			return null;
		}
		else if (proteins.size() > 1) {
			System.err.println("Too many proteins");
			return null;
		}
		else {
			Protein protein = proteins.get(0);
			for (Frame catalysis : protein.getCatalysis()) {
				boolean computationalEVCode = false;
//				boolean computationalEVCode = true; //MaizeCyc is all computational, and does not included it explicitly
				for (Object citation : catalysis.getSlotValues("CITATIONS")) {
					if (citation.toString().contains("EV-COMP")) { //We don't care how many citations, as long as one of them indicates the EC was computationally predicted
						computationalEVCode = true;
					}
				}
				if (computationalEVCode) {
					Reaction reaction = (Reaction) Reaction.load(conn, catalysis.getSlotValue("REACTION"));
					String EC = null;
					if (reaction != null) EC = reaction.getEC();
					if (EC != null && EC.length() > 0) {
						EC = EC.replace("EC-", "");//Don't allow partial EC numbers, such as 1.-.-.- or 1.2.3.-
						String p = "\\d{1,3}?.\\d{1,3}?.\\d{1,3}?.\\d{1,3}?";//EC numbers are, like IP addresses, 4 numbers of up to 3 digits each, separated by periods
						Pattern r = Pattern.compile(p);
						Matcher m = r.matcher(EC);
						if (m.find()) {
							ECList.add(EC);
						} else {
							System.err.println("partial EC? : " + EC);
						}
						
					}
				} else {
					System.err.println("Failed to find computational citation for :" + gene.getLocalID());
				}
			}
		}
		return ECList;
	}
}
