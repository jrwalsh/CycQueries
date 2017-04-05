package edu.iastate.cycqueries.queryScripts;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.iastate.javacyco.Frame;
import edu.iastate.javacyco.Gene;
import edu.iastate.javacyco.JavacycConnection;
import edu.iastate.javacyco.PtoolsErrorException;

/**
 * Written to answer a question about CornCompare paper.
 * 
 * "Would you say that the differences in your EC, pathway and reaction prediction were mostly because CornCyc and MaizeCyc use different enzyme annotation tools, 
 * or because one looked into alternative splicing products and the other did not?
 * 
 * In addition, I was wondering if you found pathway-level differences other than in photosynthesis.
 * 
 * -- Nirvana Nursimulu"
 * 
 * compareGeneEC() will grab all the gene frames (actually transcripts) from CornCyc and MaizeCyc, and only for the ones that match, print their EC assignments and whether or not
 * they are an exact match.
 * 
 * effectiveTranscripts() will determine if one "master transcript" can explain all the ECs assigned to a gene, or if 2 or more transcripts are needed to explain all
 * the EC's assigned to a gene.
 * 
 */
public class NirvanaRequest {
	
	public static void main(String[] args) {
		try {
			Long start = System.currentTimeMillis();
			
			compareGeneEC();
			effectiveTranscripts();
			
			Long stop = System.currentTimeMillis();
			Long runtime = (stop - start) / 1000;
			System.out.println("Runtime is " + runtime + " seconds.");
		}
		catch(Exception e) {
			e.printStackTrace();
			System.out.println("Caught a "+e.getClass().getName()+". Shutting down...");
		}
	}
	
	private static void compareGeneEC() throws FileNotFoundException {
		PrintWriter out = new PrintWriter("/home/jesse/Desktop/temp/out.txt");
		out.print("");
		int port = 4444;
		String organismCorn = "CORN";
		String organismMaize = "MAIZE";
		JavacycConnection conn = new JavacycConnection("localhost", port);
		
		HashMap<String, TreeSet<String>> transcriptToECMap = new HashMap<String, TreeSet<String>>();
		try {
			conn.selectOrganism(organismCorn);
			ArrayList<Frame> transcriptFrames = conn.getAllGFPInstances(Gene.GFPtype);
			for (Frame transcriptFrame : transcriptFrames) {
				String transcriptName = transcriptFrame.getCommonName();
				TreeSet<String> ECs = getGeneEC(transcriptFrame.getLocalID(), conn);
				
				//Replace any _P## with an _T## to make them consistent and therefore matchable between MaizeCyc and CornCyc
				String p = "(.*)(_P)(\\d\\d$)"; // 3 capture groups.  first group is any length of letters or numbers, followed by an _P, followed by 2 numbers and an end of line
				Pattern r = Pattern.compile(p);
				Matcher m = r.matcher(transcriptName);
				if (m.find()) {
					transcriptName = m.replaceFirst("$1_T$3");
				} else {
					//Replace any _FGP### with an _FGT### to make them consistent and therefore matchable between MaizeCyc and CornCyc
					p = "(.*)(_FGP)(\\d\\d\\d$)";
					r = Pattern.compile(p);
					m = r.matcher(transcriptName);
					if (m.find()) {
						transcriptName = m.replaceFirst("$1_FGT$3");
					}
				}
				transcriptToECMap.put(transcriptName, ECs);
			}
			
			conn.selectOrganism(organismMaize);
			transcriptFrames = conn.getAllGFPInstances(Gene.GFPtype);
			for (Frame transcriptFrame : transcriptFrames) {
				String transcriptName = transcriptFrame.getCommonName();
				TreeSet<String> ECs = getGeneEC(transcriptFrame.getLocalID(), conn);
				
				//Replace any _P## with an _T## to make them consistent and therefore matchable between MaizeCyc and CornCyc
				String p = "(.*)(_P)(\\d\\d$)"; // 3 capture groups.  first group is any length of letters or numbers, followed by an _P, followed by 2 numbers and an end of line
				Pattern r = Pattern.compile(p);
				Matcher m = r.matcher(transcriptName);
				if (m.find()) {
					transcriptName = m.replaceFirst("$1_T$3");
				} else {
					//Replace any _FGP### with an _FGT### to make them consistent and therefore matchable between MaizeCyc and CornCyc
					p = "(.*)(_FGP)(\\d\\d\\d$)";
					r = Pattern.compile(p);
					m = r.matcher(transcriptName);
					if (m.find()) {
						transcriptName = m.replaceFirst("$1_FGT$3");
					}
				}
				if (transcriptToECMap.containsKey(transcriptName)) {
					String geneName = spliceToGene(transcriptFrame.getCommonName());
					boolean match = false;
					if (transcriptToECMap.get(transcriptName).containsAll(ECs) && ECs.containsAll(transcriptToECMap.get(transcriptName))) match = true;
					if (!transcriptToECMap.get(transcriptName).isEmpty() && !ECs.isEmpty()) {
						out.println(geneName + "\t" + transcriptName + "\t" + transcriptToECMap.get(transcriptName) + "\t" + ECs + "\t" + match);
					}
				}
			}
		} catch (PtoolsErrorException e) {
			e.printStackTrace();
		}
		out.flush();
		out.close();
	}
	
	private static void effectiveTranscripts() throws FileNotFoundException {
		PrintWriter out = new PrintWriter("/home/jesse/Desktop/temp/out.txt");
		out.print("");
		int port = 4444;
		String organismCorn = "CORN";
		JavacycConnection conn = new JavacycConnection("localhost", port);
		conn.selectOrganism(organismCorn);
		
		try {
			ArrayList<Frame> transcriptFrames = conn.getAllGFPInstances(Gene.GFPtype);
			TreeMap<String, TreeSet<String>> transcriptFrameToECs = new TreeMap<String, TreeSet<String>>();
			HashMap<String, ArrayList<String>> geneToTranscriptMap = new HashMap<String, ArrayList<String>>();
			for (Frame transcriptFrame : transcriptFrames) {
				transcriptFrameToECs.put(transcriptFrame.getCommonName(), getGeneEC(transcriptFrame.getLocalID(), conn));
				String geneName = spliceToGene(transcriptFrame.getCommonName());
				if (!geneToTranscriptMap.containsKey(geneName)) geneToTranscriptMap.put(geneName, new ArrayList<String>());
				geneToTranscriptMap.get(geneName).add(transcriptFrame.getCommonName());
			}
			
			for (String geneName : geneToTranscriptMap.keySet()) {
				ArrayList<String> currentGeneECs = new ArrayList<String>();
				for (String transcriptName : geneToTranscriptMap.get(geneName)) {
					currentGeneECs.addAll(transcriptFrameToECs.get(transcriptName));
				}
				boolean hasMultipleEffectiveSpliceVariants = true;
				for (String transcriptName : geneToTranscriptMap.get(geneName)) {
					if (transcriptFrameToECs.get(transcriptName).containsAll(currentGeneECs)) hasMultipleEffectiveSpliceVariants = false; 
				}
				
				for (String transcriptName : geneToTranscriptMap.get(geneName)) {
					if (transcriptFrameToECs.get(transcriptName).containsAll(currentGeneECs)) hasMultipleEffectiveSpliceVariants = false; 
				}
				
				for (String transcriptName : geneToTranscriptMap.get(geneName)) {
					out.print(hasMultipleEffectiveSpliceVariants + "\t" + geneName + "\t" + transcriptName + "\t" + transcriptFrameToECs.get(transcriptName) + "\n");
				}
			}
			
		} catch (PtoolsErrorException e) {
			e.printStackTrace();
		}
		out.flush();
		out.close();
	}
	
	private static TreeSet<String> getGeneEC(String geneFrameID, JavacycConnection conn) {
		TreeSet<String> ecs = new TreeSet<String>();
		try {
			ArrayList<String> reactions = conn.reactionsOfGene(geneFrameID);
			for (String reaction : reactions) {
				String ec = conn.getSlotValue(reaction, "EC-Number");
				if (!ec.equalsIgnoreCase("nil") && !ec.isEmpty()) ecs.add(ec);
			}
		} catch (PtoolsErrorException e) {
			e.printStackTrace();
		}
		return ecs;
	}
	
	private static String spliceToGene(String spliceVariantName) {
		String geneName = spliceVariantName;
		geneName = geneName.replaceAll("_P..$", ""); //Remove _P## suffixs, which indicate this is a transcript
		geneName = geneName.replaceAll("_T..$", ""); //Remove _T## suffixs, which also indicate this is a transcript
		geneName = geneName.replaceAll("_FGP...$", ""); //Remove FGP#### suffixs, which also indicate this is a transcript
		geneName = geneName.replaceAll("_FGT...$", ""); //Remove FGP#### suffixs, which also indicate this is a transcript
		return geneName;
	}
}
