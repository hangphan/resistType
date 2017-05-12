package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;

import dna.Parser;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ReadWrite;
import fileIO.FileFormat;
import fileIO.TextFile;
import fileIO.TextStreamWriter;

import align2.ListNum;
import align2.Shared;
import align2.Tools;

/**
 * @author Brian Bushnell
 * @date October 8, 2014
 *
 */
public class FilterByCoverage {

	public static void main(String[] args){
		Timer t=new Timer();
		t.start();
		FilterByCoverage mb=new FilterByCoverage(args);
		mb.process(t);
	}
	
	public FilterByCoverage(String[] args){
		
		if(Parser.parseHelp(args)){
			printOptions();
			System.exit(0);
		}
		
		for(String s : args){if(s.startsWith("out=standardout") || s.startsWith("out=stdout")){outstream=System.err;}}
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");

		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		FastaReadInputStream.SPLIT_READS=false;
		stream.FastaReadInputStream.MIN_READ_LEN=1;
		Shared.READ_BUFFER_LENGTH=Tools.min(20, Shared.READ_BUFFER_LENGTH);
		Shared.capBuffers(4);
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		ReadWrite.ZIP_THREAD_DIVISOR=2;
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("null")){
				// do nothing
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				stream.FastaReadInputStream.verbose=verbose;
				ConcurrentGenericReadInputStream.verbose=verbose;
//				align2.FastaReadInputStream2.verbose=verbose;
				stream.FastqReadInputStream.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("coverage") || a.equals("cov") || a.equals("covstats") || a.equals("coverage1") || a.equals("cov1") || a.equals("covstats1")){
				covStatsAfter=b;
			}else if(a.equals("coverage0") || a.equals("cov0") || a.equals("covstats0")){
				covStatsBefore=b;
			}else if(a.equals("minc") || a.equals("mincov") || a.equals("mincoverage")){
				minCoverage=Double.parseDouble(b);
			}else if(a.equals("minp") || a.equals("minpercent")){
				minCoveredPercent=Double.parseDouble(b);
			}else if(a.equals("minr") || a.equals("minreads")){
				minReads=Long.parseLong(b);
			}else if(a.equals("minratio") || a.equals("ratio")){
				minRatio=Double.parseDouble(b);
			}else if(a.equals("minl") || a.equals("minlen") || a.equals("minlength")){
				minLength=Integer.parseInt(b);
			}else if(a.equals("appendresults") || a.equals("logappend") || a.equals("appendlog") || a.equals("appendtolog")){
				logappend=Tools.parseBoolean(b);
			}else if(a.equals("log") || a.equals("results")){
				logfile=b;
			}else if(a.equals("logheader")){
				logheader=Tools.parseBoolean(b);
			}else if(a.equals("outd") || a.equals("outdirty")){
				outdirty=b;
			}else if(parser.in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
				parser.in1=arg;
				if(arg.indexOf('#')>-1 && !new File(arg).exists()){
					parser.in1=b.replace("#", "1");
					parser.in2=b.replace("#", "2");
				}
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Download parser fields
			
			maxReads=parser.maxReads;
			
			overwrite=parser.overwrite;
			append=parser.append;
			
			if(parser.minReadLength>0){minLength=parser.minReadLength;}
			
			in1=parser.in1;
			qfin1=parser.qfin1;

			outclean=parser.out1;
			qfoutclean=parser.qfout1;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		name=ReadWrite.stripToCore(in1);
		
		if(!ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(outclean!=null && outclean.equalsIgnoreCase("null")){outclean=null;}
		if(outdirty!=null && outdirty.equalsIgnoreCase("null")){outdirty=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, outclean, outdirty)){
			outstream.println((outclean==null)+", "+outclean+", "+(outdirty==null)+", "+outdirty);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+outclean+", "+outdirty+"\n");
		}
		
		{
			byte qin=Parser.qin, qout=Parser.qout;
			if(qin!=-1 && qout!=-1){
				FASTQ.ASCII_OFFSET=qin;
				FASTQ.ASCII_OFFSET_OUT=qout;
				FASTQ.DETECT_QUALITY=false;
			}else if(qin!=-1){
				FASTQ.ASCII_OFFSET=qin;
				FASTQ.DETECT_QUALITY=false;
			}else if(qout!=-1){
				FASTQ.ASCII_OFFSET_OUT=qout;
				FASTQ.DETECT_QUALITY_OUT=false;
			}
		}

		ffoutclean=FileFormat.testOutput(outclean, FileFormat.FASTA, extout, true, overwrite, append, false);
		ffoutdirty=FileFormat.testOutput(outdirty, FileFormat.FASTA, extout, true, overwrite, append, false);

		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, extin, true, true);
		ffCov0=FileFormat.testInput(covStatsBefore, FileFormat.TEXT, ".txt", true, false);
		ffCov1=FileFormat.testInput(covStatsAfter, FileFormat.TEXT, ".txt", true, false);
	}
	
	void process(Timer t){

		final HashMap<String, CovStatsLine> cslMap0=new HashMap<String, CovStatsLine>(1024);
		final HashMap<String, CovStatsLine> cslMap1=new HashMap<String, CovStatsLine>(1024);
		if(ffCov0!=null){
			TextFile tf=new TextFile(ffCov0);
			int i=0;
			for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
				if(!s.startsWith("#") && i>0){
					CovStatsLine csl=new CovStatsLine(s);
					CovStatsLine old=cslMap0.put(csl.id, csl);
					assert(old==null);
				}
				i++;
			}
			tf.close();
		}
		if(ffCov1!=null){
			TextFile tf=new TextFile(ffCov1);
			int i=0;
			for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
				if(!s.startsWith("#") && i>0){
					CovStatsLine csl=new CovStatsLine(s);
					CovStatsLine old=cslMap1.put(csl.id, csl);
					assert(old==null);
				}
				i++;
			}
			tf.close();
		}
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, false, false, ffin1, null, qfin1, null);
			if(verbose){outstream.println("Started cris");}
			cris.start(); //4567
		}
		assert(!cris.paired());
		
		final ConcurrentReadOutputStream rosClean;
		if(outclean!=null){
			final int buff=4;
			
			assert(!outclean.equalsIgnoreCase(in1) && !outclean.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			rosClean=ConcurrentReadOutputStream.getStream(ffoutclean, null, qfoutclean, null, buff, null, false);
			rosClean.start();
		}else{rosClean=null;}
		
		final ConcurrentReadOutputStream rosDirty;
		if(outdirty!=null){
			final int buff=4;
			
			assert(!outdirty.equalsIgnoreCase(in1) && !outdirty.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			rosDirty=ConcurrentReadOutputStream.getStream(ffoutdirty, null, qfoutdirty, null, buff, null, false);
			rosDirty.start();
		}else{rosDirty=null;}
		
		long readsProcessed=0;
		long basesProcessed=0;
		
		long readsOut=0;
		long basesOut=0;
		
		long readsFiltered=0;
		long basesFiltered=0;
		
		final TextStreamWriter tsw=(logfile==null ? null : new TextStreamWriter(logfile, (overwrite && !logappend), logappend, true));
//		System.err.println("***** overwrite="+overwrite+", logappend="+logappend+", combined="+(overwrite && !logappend));
		if(tsw!=null){
			tsw.start();
			if(logheader){tsw.print("#assembly\tcontig\tcontam\tlength\tavgFold\treads\tpercentCovered"+(ffCov0==null ? "" : "\tavgFold0\treads0\tnormRatio")+"\n");}
		}
		
		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
//			outstream.println("Fetched "+reads);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}

			while(reads!=null && reads.size()>0){

				final ArrayList<Read> cleanList=new ArrayList<Read>(reads.size());
				final ArrayList<Read> dirtyList=new ArrayList<Read>(reads.size());
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					assert(r1.mate==null);
					
					final int initialLength1=r1.length();
					
					readsProcessed++;
					basesProcessed+=initialLength1;
					
					final double ratio;
					final boolean contam;
					
					final CovStatsLine csl0=cslMap0.get(r1.id);
					final CovStatsLine csl1=cslMap1.get(r1.id);
					if(csl1!=null){
						ratio=(csl0==null ? 0 : csl0.avgFold/Tools.max(0.01, csl1.avgFold));

						if(csl1.reads()<minReads || csl1.length<minLength || csl1.coveredPercent()<minCoveredPercent){
							contam=true;
						}else if((csl1.avgFold<minCoverage && ratio>minRatio) || csl1.avgFold<0.5){
							contam=true;
						}else{
							contam=false;
						}
					}else{
						contam=true;
						ratio=0;
					}
					
					
					if(!contam){
						cleanList.add(r1);
						readsOut++;
						basesOut+=initialLength1;
					}else{
						dirtyList.add(r1);
						readsFiltered++;
						basesFiltered+=initialLength1;
					}
					if(tsw!=null && (r1.length()>=minLength || PRINT_SHORT_CONTIG_RESULTS)){
						if(csl1==null){
							if(ffCov0==null){
								tsw.print(String.format("%s\t%s\t%s\t%d\t%.2f\t%d\t%.2f\n", name, r1.id, contam ? "1" : "0", r1.length(), 0, 0, 0));
							}else{
								tsw.print(String.format("%s\t%s\t%s\t%d\t%.2f\t%d\t%.2f\t%.2f\t%d\t%.2f\n", 
										name, r1.id, contam ? "1" : "0", r1.length(), 0, 0, 0, 0, 0, 0));
							}
							
						}else if(csl0==null){
							tsw.print(String.format("%s\t%s\t%s\t%d\t%.2f\t%d\t%.2f\n", name, csl1.id, contam ? "1" : "0", csl1.length,
									csl1.avgFold, csl1.plusReads+csl1.minusReads, csl1.coveredPercent()));
						}else{
							tsw.print(String.format("%s\t%s\t%s\t%d\t%.2f\t%d\t%.2f\t%.2f\t%d\t%.2f\n", name, csl1.id, contam ? "1" : "0", csl1.length,
									csl1.avgFold, csl1.plusReads+csl1.minusReads, csl1.coveredPercent(), csl0.avgFold, csl0.plusReads+csl0.minusReads, ratio));
						}
					}
				}

				if(rosClean!=null){rosClean.add(cleanList, ln.id);}
				if(rosDirty!=null){rosDirty.add(dirtyList, ln.id);}

				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
		}
		
		errorState|=ReadWrite.closeStreams(cris, rosClean, rosDirty);
		if(tsw!=null){errorState|=tsw.poisonAndWait();}
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);
		
		outstream.println("Time:               "+t);
		outstream.println("Reads In:           "+readsProcessed+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases In:           "+basesProcessed+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		outstream.println("Reads Out:          "+readsOut);
		outstream.println("Bases Out:          "+basesOut);
		outstream.println("Reads Filtered:     "+readsFiltered);
		outstream.println("Bases Filtered:     "+basesFiltered);
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/*--------------------------------------------------------------*/
	
	private void printOptions(){
		assert(false) : "printOptions: TODO";
//		outstream.println("Syntax:\n");
//		outstream.println("java -ea -Xmx512m -cp <path> jgi.ReformatReads in=<infile> in2=<infile2> out=<outfile> out2=<outfile2>");
//		outstream.println("\nin2 and out2 are optional.  \nIf input is paired and there is only one output file, it will be written interleaved.\n");
//		outstream.println("Other parameters and their defaults:\n");
//		outstream.println("overwrite=false  \tOverwrites files that already exist");
//		outstream.println("ziplevel=4       \tSet compression level, 1 (low) to 9 (max)");
//		outstream.println("interleaved=false\tDetermines whether input file is considered interleaved");
//		outstream.println("fastawrap=80     \tLength of lines in fasta output");
//		outstream.println("qin=auto         \tASCII offset for input quality.  May be set to 33 (Sanger), 64 (Illumina), or auto");
//		outstream.println("qout=auto        \tASCII offset for output quality.  May be set to 33 (Sanger), 64 (Illumina), or auto (meaning same as input)");
//		outstream.println("outsingle=<file> \t(outs) Write singleton reads here, when conditionally discarding reads from pairs.");
	}
	
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String covStatsBefore=null;
	private String covStatsAfter=null;
	private String name=null;
	
	private String qfin1=null;

	private String outclean=null;
	private String outdirty=null;

	private String qfoutclean=null;
	private String qfoutdirty=null;
	
	private String extin=null;
	private String extout=null;
	
	/*--------------------------------------------------------------*/
	
	private long maxReads=-1;

	/** Scaffolds shorter than this will be discarded. */
	private int minLength=0;
	/** Scaffolds with fewer mapped reads will be discarded. */
	private long minReads=0;
	/** Scaffolds with lower average coverage will be discarded. */
	private double minCoverage=5;
	/** Scaffolds with a lower percent of covered bases will be discarded. */
	private double minCoveredPercent=40;
	/** Scaffolds will NOT be discarded based on low coverage unless the coverage dropped by at least this factor. */
	private double minRatio=0;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffCov0;
	private final FileFormat ffCov1;

	private final FileFormat ffoutclean;
	private final FileFormat ffoutdirty;
	
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	private boolean logappend=false;
	private String logfile=null;
	private boolean logheader=true;
	private static boolean PRINT_SHORT_CONTIG_RESULTS=false;
	
}