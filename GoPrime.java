package goprime;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;

public class GoPrime {

    public static String primers[],  seqs[], compSeqs[], seqNames[], pNames[];
    public static String pFilename, sFilename, outFilename, head, header, tTx;
    
    public static int scores[][][], pLens[], mLen, count;
    
    public static double minPm, minPmI, minPb;
    
    public static char bases[][], revCodes[][], allow[], amb[][], comps[][], compsBase[], ts[][], muts[][];
    
    public static double primerCts[], primerCts2[], probeCts[], intraPrimer, interPrimer, intraProbe, interProbe, deltaLod, max14nt;
    
    public static DecimalFormat df;
    
    public static boolean fullOut, seqOut, success;
    
    public static int bestProbe, bestProbeSt, bestProbeEn;
    public static double bestProbePer;
    
    public static String noHits=""+System.getProperty("line.separator");
    
    
    public static void main(String[] args) {
        
        System.out.println("GoPrime Started"+System.getProperty("line.separator"));

        df = new DecimalFormat();
        df.setMaximumFractionDigits(2);
        
        if(args.length==2) {
            pFilename=args[0];
            sFilename=args[1];
        }
        else {
            System.out.println("Error - incorrect usage: "+args.length +" parameters were given");
            System.out.println("Correct usage = java -jar GoPrime.jar primerSeqFilename.fasta targetSeqFilename.fasta");
            System.out.println("Exiting..."+System.getProperty("line.separator"));
            System.exit(0);           
        }

        System.out.println("Primer Sequence File = "+pFilename);
        System.out.println("Target Sequence File = "+sFilename);
        
        //PRINT Output file names
        
        setUpBases();
        setUpCts();
        
        primers=new String[6];
        for(int i=0;i<primers.length;i++) {
            primers[i]="";
        }
        
        pNames=new String[primers.length];
        pNames[0]="5'-Fwd-3'-";//matches onto 3'-5' strand
        pNames[1]="5'-Probe-3'-";//matches onto 3'-5' strand
        pNames[2]="5'-Rev-3'-";//matches onto 3'-5' strand
        pNames[3]="3'-Fwd-5'-";//matches onto 5'-3' strand
        pNames[4]="3'-Probe-5'-";//matches onto 5'-3' strand
        pNames[5]="3'-Rev-5'-";//matches onto 5'-3' strand
        
        //Primers must be FASTA, 5'-3' format
        File inFile = new File(pFilename);
        count=0;
        try {
            BufferedReader input =  new BufferedReader(new FileReader(inFile));

            try {
                String line = null;

                while (( line = input.readLine()) != null) {
                    if(line.indexOf(">")==0) {
                        count++;
                        
                        //Primer names
                        pNames[count-1]+=line.substring(1, line.length());
                        pNames[count-1+3]+=line.substring(1, line.length());
                        
                        if(count>3) {
                            System.out.println("Error - more than 3 sequences in primer file");
                            System.out.println("Exiting...");
                            System.exit(0);
                        }
                    }
                    else {
                        primers[count-1]+=line;
                    }
                }
            }
            finally {
                input.close();
            }
        }
        catch (IOException ex) {
            ex.printStackTrace();
        }
        
        if(count==2) {
            System.out.println("Only two sequences in primers file - creating blank probe sequence");
            primers[2]=primers[1];
            primers[5]=primers[1];
            primers[1]="NNNNNNNNNNNNNNNN";
            primers[4]=primers[1];
        }
        else if(count!=3) {
            System.out.println("Error - abnormal number of sequences in primer file (expected 2 primers and/or 1 probe");
            System.out.println("Exiting...");
            System.exit(0);
        }
        
        System.out.println("Forward Primer Seq = "+primers[0]+" length = "+primers[0].length());
        System.out.println("Probe Seq = "+primers[1]+" length = "+primers[1].length());
        System.out.println("Reverse Primer Seq = "+primers[2]+" length = "+primers[2].length());
        
        //CheckPrimers
        for(int p=0;p<3;p++) {
            //Technically GAPs are in the 'allowables' and get removed but this is for future development
            if(primers[p].indexOf("-")>=0) {
                System.out.println("Error - primer sequence has a gap in it = "+primers[p]);
                System.out.println("Exiting...");
                System.exit(0);
            }
            
            primers[p]=checkSeqs(pNames[p],primers[p]);
        }
        
        //ReversePrimers - NOT reverse-complement
        //The original 5'3' primers are checked against the complement target sequence which are in 3'-5'
        //Then check the reverse primers (now 3'-5') against the original 5'-3' target seqs
        for(int p=0;p<3;p++) {
            primers[p+3]=reverse(primers[p]);
        }
             
        pLens=new int[primers.length];
        for(int p=0;p<primers.length;p++) {
            pLens[p]=primers[p].length();
        }
        
        //Combined primer length - needed for combined mismatch freq
        mLen=primers[0].length()+primers[2].length();
        
        inFile = new File(sFilename);
        count=0;
        try {
            BufferedReader input =  new BufferedReader(new FileReader(inFile));

            try {
                String line = null;

                while (( line = input.readLine()) != null) {
                    if(line.indexOf(">")==0) {
                        count++;
                    }
                }
            }
            finally {
                input.close();
            }
        }
        catch (IOException ex) {
            ex.printStackTrace();
        }
        
        if(count==0) {
            System.out.println("Error - no seqs found - check file is FASTA format?");
            System.out.println("Exiting...");
            System.exit(0);
        }
        
        seqs=new String[count];
        seqNames=new String[count];
        for(int i=0;i<seqs.length;i++) {
            seqs[i]="";
            seqNames[i]="";
        }
        
        count=0;
        try {
            BufferedReader input =  new BufferedReader(new FileReader(inFile));

            try {
                String line = null;

                while (( line = input.readLine()) != null) {
                    if(line.indexOf(">")==0) {
                        count++;
                        seqNames[count-1]=line;
                    }
                    else {
                        seqs[count-1]+=line;
                    }
                }
            }
            finally {
                input.close();
            }
        }
        catch (IOException ex) {
            ex.printStackTrace();
        }
        
        System.out.println("");
        System.out.println(count+" sequences found in target sequence file");
        
        for(int i=0;i<seqs.length;i++) {
            seqs[i]=checkSeqs(seqNames[i], seqs[i]);
        }

        compSeqs=new String[seqs.length];
        for(int i=0;i<seqs.length;i++) {
            compSeqs[i]=complement(seqs[i]);
        }
        
        System.out.println("");
        
        try {
            outFilename=sFilename+"_out.txt";

            FileWriter fstreamOut = new FileWriter(outFilename);
            BufferedWriter outOut = new BufferedWriter(fstreamOut);
            
            outFilename=sFilename+"_cts.txt";

            FileWriter fstreamCT = new FileWriter(outFilename);
            BufferedWriter outCT = new BufferedWriter(fstreamCT);
        
            //Loop through each sequence in turn
            for(int s=0;s<compSeqs.length;s++) {
                success=false;
                
                try {
                    //outFilename=sFilename.substring(0, sFilename.lastIndexOf("."))+"_"+s+"_matrix.txt";
                    outFilename=sFilename+"_"+(s+1)+"_matrix.txt";

                    tTx="Evaluating seq "+(s+1)+" "+seqNames[s];
                    System.out.println(tTx);
                    outOut.write(tTx+System.getProperty("line.separator"));

                    //Clear results - 6 for F/P/R*2, 14 for all the different mutation counts - plus one for Ns now
                    scores=new int[compSeqs[s].length()][6][16];

                    //F/P/R primers in ForwardDirection checked against the ComplementSeqs
                    for(int p=0;p<3;p++) {
                        for(int i=0;i<=compSeqs[s].length()-pLens[p];i++) {

                            //this is where to record the result - want it on the 3' end of the primer binding site
                            int pos=i+pLens[p]-1;

                            for(int j=0;j<pLens[p];j++) {

                                if(compSeqs[s].charAt(i+j)=='N') 
                                    scores[pos][p][14]++;//number of Ns
                                
                                int mut=checkMutation(primers[p].charAt(j),compSeqs[s].charAt(i+j));

                                if(mut>0) {
                                    scores[pos][p][0]++;//total mismatches

                                    if(mut==1) {
                                        scores[pos][p][1]++;//total ts
                                    }
                                    else if(mut==2) {
                                        scores[pos][p][2]++;//total tv
                                    }

                                    if(j==pLens[p]-1) {
                                        if(mut==1)
                                            scores[pos][p][3]++;//ts at nt1
                                        else if(mut==2)
                                            scores[pos][p][4]++;//tv at nt1

                                        scores[pos][p][11]++;//mismatches nt1-4
                                    }

                                    else if(j==pLens[p]-2) {
                                        if(mut==1)
                                            scores[pos][p][5]++;//ts at nt2
                                        else if(mut==2)
                                            scores[pos][p][6]++;//tv at nt2

                                        scores[pos][p][11]++;//mismatches nt1-4
                                    }

                                    else if(j==pLens[p]-3) {
                                        if(mut==1)
                                            scores[pos][p][7]++;//ts at nt3
                                        else if(mut==2)
                                            scores[pos][p][8]++;//tv at nt3

                                        scores[pos][p][11]++;//mismatches nt1-4
                                    }

                                    else if(j==pLens[p]-4) {
                                        if(mut==1)
                                            scores[pos][p][9]++;//ts at nt4
                                        else if(mut==2)
                                            scores[pos][p][10]++;//tv at nt4

                                        scores[pos][p][11]++;//mismatches nt1-4
                                    }
                                    
                                    if(j<4) {
                                        scores[pos][p][15]++;//mismatches at 5' 1-4nt
                                    }
                                        
                                }//end of mut>0
                            }//end or primer characters loop j
                        }//end of sequence characters loop i
                    }//end of F/P/R primer (forward) loop p

                    //F/P/R primers in ReverseDirection checked against Seqs
                    for(int p=3;p<6;p++) {
                        for(int i=0;i<=seqs[s].length()-pLens[p];i++) {

                            int pos=i;

                            for(int j=0;j<pLens[p];j++) {

                                if(seqs[s].charAt(i+j)=='N') 
                                    scores[pos][p][14]++;
                                
                                int mut=checkMutation(primers[p].charAt(j),seqs[s].charAt(i+j));

                                if(mut>0) {
                                    scores[pos][p][0]++;//total mismatches

                                    if(mut==1) {
                                        scores[pos][p][1]++;//total ts
                                    }
                                    else if(mut==2) {
                                        scores[pos][p][2]++;//total tv
                                    }

                                    if(j==0) {
                                        if(mut==1)
                                            scores[pos][p][3]++;//ts at nt1
                                        else if(mut==2)
                                            scores[pos][p][4]++;//tv at nt1

                                        scores[pos][p][11]++;//mismatches nt1-4
                                    }

                                    else if(j==1) {
                                        if(mut==1)
                                            scores[pos][p][5]++;//ts at nt2
                                        else if(mut==2)
                                            scores[pos][p][6]++;//tv at nt2

                                        scores[pos][p][11]++;//mismatches nt1-4
                                    }

                                    else if(j==2) {
                                        if(mut==1)
                                            scores[pos][p][7]++;//ts at nt3
                                        else if(mut==2)
                                            scores[pos][p][8]++;//tv at nt3

                                        scores[pos][p][11]++;//mismatches nt1-4
                                    }

                                    else if(j==3) {
                                        if(mut==1)
                                            scores[pos][p][9]++;//ts at nt4
                                        else if(mut==2)
                                            scores[pos][p][10]++;//tv at nt4

                                        scores[pos][p][11]++;//mismatches nt1-4
                                    }
                                    
                                    if(j>=pLens[p]-4){
                                        scores[pos][p][15]++;//mismactehs at 5' 1-4nt
                                    }
                                }//end of mut>0
                            }//end or primer characters loop j
                        }//end of sequence characters loop i
                    }//end of F/P/R primer loop p

                    //Set the overhangs at the ends to max mismatches - only records score at the 3' position of primer so empty values at the ends
                    for(int i=0;i<scores.length;i++) {//i=position
                        for(int j=0;j<scores[i].length;j++) {//j=primer
                            if(j<3) {
                                //Matching Forward Direction Primers to 3'-5' stand so empty vals at left 3' end
                                if(i<pLens[j]-1) {
                                    scores[i][j][0]=pLens[j];
                                }
                            }
                            else {
                                //Matching Reverse Direction Primers to 5'-3' stand so empty vals at right 3' end
                                if(i>compSeqs[s].length()-pLens[j]) {
                                    scores[i][j][0]=pLens[j];
                                }
                            }
                        }
                    }

                    //Rule out positions based on % mismatch and number of mismatches at 4nt 3' end
                    for(int i=0;i<scores.length;i++) {
                        for(int j=0;j<scores[i].length;j++) {
                            scores[i][j][13]=0;//0=pass, >0 is fail, set to pass initially

                            //check Ns - bit curde - what to discount regions with high N
                            double ns=(double)scores[i][j][14]/(double)pLens[j];
                            if(ns>=minPb | ns>=minPm)//
                                scores[i][j][13]+=4;
                            
                            //probe
                            if(j==1 | j==4) {
                                double perc=(double)scores[i][j][0]/(double)primers[j].length();
                                if(perc>=minPb) {
                                   scores[i][j][12]=1;//flag for failed % test 
                                   scores[i][j][13]+=2; 
                                }
                            }
                            //primer
                            else {
                                //if more than 2 mismatches in 1-4nt at 3'
                                if(scores[i][j][11]>max14nt) {
                                   scores[i][j][13]+=1; 
                                }
                                //use mLen (combined F and R length) initially to find initial candidates - filter later
                                double perc=(double)scores[i][j][0]/(double)mLen;
                                double percI=(double)scores[i][j][0]/(double)pLens[j];
                                
                                if(perc>=minPm | percI>=minPmI) {
                                   scores[i][j][12]=1;//flag for failed % test 
                                   scores[i][j][13]+=2; 
                                }
                            }
                        }
                    }//end of finding positions that prime loop

                    //count candidate positions
                    int can[]=new int[6];

                    for(int i=0;i<scores[0].length;i++) {
                        count=0;

                        for(int j=0;j<scores.length;j++) {
                            if(scores[j][i][13]==0) {
                                double mis=100-(double)scores[j][i][0]/(double)pLens[i]*100;
                                double mis2=mis;
                                if(i!=1 & i!=4)
                                    mis2=100-(double)scores[j][i][0]/(double)(mLen)*100;

                                tTx=(j+1)+" position is a candidate for "+pNames[i]+" "+df.format(mis)+"% "+df.format(mis2)+"% [%Match %MatchPair], "+scores[j][i][0]+" "+scores[j][i][11]+" [TotMis Tot1-4]";
                                System.out.println(tTx);
                                outOut.write(tTx+System.getProperty("line.separator"));
                                count++;
                            }
                        }
                        can[i]=count;
                    }

                    //RJO - at this point the probe is not checked if it is inbetween the two 
                    //Technically probes could go in any set/direction as long as between the F/R primers
                    tTx=can[0]+"-["+(can[1]+can[4])+"]-"+can[5]+" Fwd-[Probe]-Rev individual candidate primer/probe positions found in expected orientation";
                    System.out.println(tTx);
                    outOut.write(tTx+System.getProperty("line.separator"));
                    
                    tTx=can[3]+"-["+(can[1]+can[4])+"]-"+can[2]+" Fwd-[Probe]-Rev individual candidate primer/probe positions found in opposite orientation";
                    System.out.println(tTx);
                    outOut.write(tTx+System.getProperty("line.separator"));
                    
                    int oCount=0;

                    int fSt=0, fEn=0, pSt=0, pEn=0, rSt=0, rEn=0;

                    //scores only stores data related to this seq (s) - so i is length of seq s
                    for(int i=0;i<scores.length;i++) {

                        fSt=fEn=pSt=pEn=rSt=rEn=0;

                        //Forward(5'-3') of Forward primer is good
                        if(scores[i][0][13]==0) { 

                            fSt=i-pLens[0]+1;
                            fEn=i;

                            for(int j=i+1;j<scores.length;j++) {//only plus one (rather than primer length) as rev primer is in reverse direction

                                //Reverse (3'-5') of Reverse primer is good
                                if(scores[j][5][13]==0) {

                                    bestProbe=1000000;
                                    bestProbeSt=bestProbeEn=0;
                                    bestProbePer=0;
                                    
                                    rSt=j+pLens[5]-1;//chnaged from [2] to [5] - should be same
                                    rEn=j;

                                    String fMes="";
                                    boolean oPass=true;

                                    //Check combined mismatch and % OK
                                    int mis14=scores[i][0][11]+scores[j][5][11];//[0]is FwdPrimer in 5'-3' and [5] is RevPrimer in 3'-5'
                                    if(mis14>max14nt) {
                                        oPass=false;
                                        fMes+="(Rule1) Failed "+mis14+">"+max14nt+" mismatches in 3' 1-4nts of the pair ";
                                    }
                                    double pairM=(double)(scores[i][0][0]+scores[j][5][0])/(double)mLen;
                                    double pairM2=100-pairM*100;
                                    if(pairM>=minPm) {
                                        oPass=false;
                                        fMes+="(Rule2) Failed "+df.format(pairM2)+"%<"+(minPm*100)+"% matches of the pair ";
                                    } 

                                    boolean test=false;

                                    //Forward (5'-3') of the Probe
                                    //for(int k=i+pLens[1];k<j;k++) {//old way - with no probe overlaps
                                    for(int k=i+pLens[1];k<=rSt;k++) {//rSt is the end really as rev primer is 3'-5'
                                       
                                        if(scores[k][1][0]<bestProbe) {
                                            bestProbe=scores[k][1][0];
                                            bestProbeSt=k-pLens[1]+1;
                                            bestProbeEn=k;
                                            bestProbePer=100-(double)scores[k][1][0]/(double)pLens[1]*100;
                                        }
                                        //if fwd probe is good
                                        if(scores[k][1][13]==0) {

                                            pSt=k-pLens[1]+1;
                                            pEn=k;

                                            if(oPass) {
                                                oCount++;
                                                //System.out.println("S"+oCount+" RT-PCR Success Fwd="+(i-pLens[0]+1)+"/"+i+" Probe="+(k-pLens[1]+1)+"/"+k+" Rev="+(j+pLens[5]-1)+"/"+j);

                                                double cts[]=calculateCT(fEn,0,pEn,1,rEn,5);
                                                tTx="Set-"+oCount+" RT-PCR Success Fwd="+(fSt+1)+"/"+(fEn+1)+" Probe="+(pSt+1)+"/"+(pEn+1)+" Rev="+(rSt+1)+"/"+(rEn+1);
                                                System.out.println(tTx);
                                                outOut.write(tTx+System.getProperty("line.separator"));
                                                
                                                tTx="deltaCT="+df.format(cts[0])+" [fCT="+df.format(cts[1])+" pCT="+df.format(cts[2])+" rCT="+df.format(cts[3])+"]";
                                                System.out.println(tTx);
                                                outOut.write(tTx+System.getProperty("line.separator"));
                                                
                                                outCT.write(seqNames[s]+"\t"+df.format(cts[0])+"\t"+df.format(cts[4])+"\t"+"Fwd[CT="+df.format(cts[1])+"]="+(fSt+1)+"/"+(fEn+1)+"\tProbe=[CT="+df.format(cts[2])+"]="+(pSt+1)+"/"+(pEn+1)+"\tRev=[CT="+df.format(cts[3])+"]="+(rSt+1)+"/"+(rEn+1)+"\t"+(rSt-fSt)+System.getProperty("line.separator"));
                                                success=true;
                                            }
                                            else {
                                                oCount++;
                                                //System.out.println("S"+oCount+" Failed Fwd="+(i-pLens[0]+1)+"/"+i+" Probe="+(k-pLens[1]+1)+"/"+k+" Rev="+(j+pLens[5]-1)+"/"+j+" "+fMes);
         
                                                tTx="Set-"+oCount+" Failed Fwd="+(fSt+1)+"/"+(fEn+1)+" Probe="+(pSt+1)+"/"+(pEn+1)+" Rev="+(rSt+1)+"/"+(rEn+1)+" "+fMes;
                                                System.out.println(tTx);
                                                outOut.write(tTx+System.getProperty("line.separator"));
                                            
                                            }
                                            test=true;
                                        }
                                    }
                       

                                    //j is rEn
                                    // minus the length of probe gives 
                                    
                                    //Reverse (3'-5') of the Probe
                                    //for(int k=i+1;k<j-pLens[1];k++) {//old way - with no probe overlaps
                                    for(int k=fSt;k<j-pLens[1];k++) {    
                                        
                                        if(scores[k][4][0]<bestProbe) {
                                            bestProbe=scores[k][4][0];
                                            bestProbeSt=k+pLens[4]+1;
                                            bestProbeEn=k;
                                            bestProbePer=100-(double)scores[k][4][0]/(double)pLens[4]*100;
                                        }
                                        
                                        if(scores[k][4][13]==0) {

                                            pSt=k+pLens[4]+1;
                                            pEn=k;

                                            if(oPass) {
                                                oCount++;
                                                //System.out.println("S"+oCount+" RT-PCR Success Fwd="+(i-pLens[0]+1)+"/"+i+" Probe="+k+"/"+(k-pLens[4]+1)+" Rev="+(j+pLens[5]-1)+"/"+j);

                                                double cts[]=calculateCT(fEn,0,pEn,4,rEn,5);
                                                tTx="Set-"+oCount+" RT-PCR Success Fwd="+(fSt+1)+"/"+(fEn+1)+" Probe="+(pSt+1)+"/"+(pEn+1)+" Rev="+(rSt+1)+"/"+(rEn+1);
                                                System.out.println(tTx);
                                                outOut.write(tTx+System.getProperty("line.separator"));
                                                
                                                tTx="deltaCT="+df.format(cts[0])+" [fCT="+df.format(cts[1])+" pCT="+df.format(cts[2])+" rCT="+df.format(cts[3])+"]";
                                                System.out.println(tTx);
                                                outOut.write(tTx+System.getProperty("line.separator"));
                                                
                                                outCT.write(seqNames[s]+"\t"+df.format(cts[0])+"\t"+df.format(cts[4])+"\t"+"Fwd[CT="+df.format(cts[1])+"]="+(fSt+1)+"/"+(fEn+1)+"\tProbe=[CT="+df.format(cts[2])+"]="+(pSt+1)+"/"+(pEn+1)+"\tRev=[CT="+df.format(cts[3])+"]="+(rSt+1)+"/"+(rEn+1)+"\t"+(rSt-fSt)+System.getProperty("line.separator"));
                                                success=true;
                                            }
                                            else {
                                                oCount++;
                                                //System.out.println("S"+oCount+" Failed Fwd="+(i-pLens[0]+1)+"/"+i+" Probe="+k+"/"+(k-pLens[4]+1)+" Rev="+(j+pLens[5]-1)+"/"+j+" "+fMes);                                     
                                                tTx="Set-"+oCount+" Failed Fwd="+(fSt+1)+"/"+(fEn+1)+" Probe="+(pSt+1)+"/"+(pEn+1)+" Rev="+(rSt+1)+"/"+(rEn+1)+" "+fMes;
                                                System.out.println(tTx);
                                                outOut.write(tTx+System.getProperty("line.separator"));
                                            }
                                            test=true;
                                        }
                                    }

                                    //if didn't find a passed probe site
                                    if(!test) {
                                        oPass=false;

                                        fMes+="(Rule 3) Failed as without [effective] probe ";

                                        if(!fMes.equals("")) {
                                            oCount++;
                                            //System.out.println("S"+oCount+" Failed Fw="+(i-pLens[0]+1)+"/"+i+" Rv="+(j+pLens[2]-1)+"/"+j+" "+fMes);
                                            tTx="Set-"+oCount+" Failed Fw="+(fSt+1)+"/"+(fEn+1)+" Rv="+(rSt+1)+"/"+(rEn+1)+" BestFailedProbe="+df.format(bestProbePer)+"% "+(bestProbeSt+1)+"/"+(bestProbeEn+1)+" "+fMes;
                                            System.out.println(tTx);
                                            outOut.write(tTx+System.getProperty("line.separator"));
                                        }
                                    }
                                }//rev primer good
                            }//(rev) rev primer loop
                        }//fwd primer good
                    }//(fwd) fwd primer loop

                    //Now we look for the opposite sets
                    for(int i=0;i<scores.length;i++) {
                        //Forward (5'-3') of the Reverse primer is good
                        if(scores[i][2][13]==0) { 

                            rSt=i-pLens[2]+1;
                            rEn=i;

                            for(int j=i+1;j<scores.length;j++) {//only plus one as fwd primer is in reverse direction

                                //Reverse (3'-5')of the Forward primer is good
                                if(scores[j][3][13]==0) {

                                    bestProbe=1000000;
                                    bestProbeSt=bestProbeEn=0;
                                    bestProbePer=0;
                                    
                                    fSt=j+pLens[3]-1;
                                    fEn=j;

                                    String fMes="";
                                    boolean oPass=true;
                                    
                                    //Check combined mismatch and % OK
                                    int mis14=scores[i][2][11]+scores[j][3][11];
                                    if(mis14>max14nt) {
                                        oPass=false;
                                        fMes+="(Rule1) Failed "+mis14+">"+max14nt+" mismatches in 3' 1-4nts of the pair ";
                                    }
                                    double pairM=(double)(scores[i][2][0]+scores[j][3][0])/(double)mLen;
                                    double pairM2=100-pairM*100;
                                    if(pairM>=minPm) {
                                        oPass=false;
                                        fMes+="(Rule2) Failed "+df.format(pairM2)+"%<82.05% matches of the pair ";
                                    } 
                               
                                    boolean test=false;
                                    //Forward (5'-3') of the Probe
                                    //for(int k=i+pLens[1];k<j;k++) {//old way - with no probe overlaps
                                    for(int k=i+pLens[1];k<=fSt;k++) {
                                        
                                        if(scores[k][1][0]<bestProbe) {
                                            bestProbe=scores[k][1][0];
                                            bestProbeSt=k-pLens[1]+1;
                                            bestProbeEn=k;
                                            bestProbePer=100-(double)scores[k][1][0]/(double)pLens[1]*100;
                                        }
                                        
                                        pSt=k-pLens[1]+1;
                                        pEn=k;

                                        if(scores[k][1][13]==0) {
                                            if(oPass) {
                                                oCount++;
                                                //System.out.println("S"+oCount+" RT-PCR Success Rev="+(i-pLens[2]+1)+"/"+i+" Probe="+k+"/"+(k-pLens[1]+1)+" Fwd="+(j+pLens[3]-1)+"/"+j);

                                                double cts[]=calculateCT(fEn,3,pEn,1,rEn,2);
                                                tTx="oSet-"+oCount+" RT-PCR Success Rev="+(rSt+1)+"/"+(rEn+1)+" Probe="+(pSt+1)+"/"+(pEn+1)+" Fwd="+(fSt+1)+"/"+(fEn+1);
                                                System.out.println(tTx);
                                                outOut.write(tTx+System.getProperty("line.separator"));
                                                
                                                tTx="deltaCT="+df.format(cts[0])+" [fCT="+df.format(cts[1])+" pCT="+df.format(cts[2])+" rCT="+df.format(cts[3])+"]";
                                                System.out.println(tTx);
                                                outOut.write(tTx+System.getProperty("line.separator"));
                                                
                                                outCT.write(seqNames[s]+"\t"+df.format(cts[0])+"\t"+df.format(cts[4])+"\t"+"Fwd[CT="+df.format(cts[1])+"]="+(fSt+1)+"/"+(fEn+1)+"\tProbe=[CT="+df.format(cts[2])+"]="+(pSt+1)+"/"+(pEn+1)+"\tRev=[CT="+df.format(cts[3])+"]="+(rSt+1)+"/"+(rEn+1)+"\t"+(fSt-rSt)+System.getProperty("line.separator"));
                                                success=true;
                                            }
                                            else {
                                                oCount++;
                                                //System.out.println("S"+oCount+" Failed Rev="+(i-pLens[2]+1)+"/"+i+" Probe="+k+"/"+(k-pLens[0]+1)+" Fwd="+(j+pLens[3]-1)+"/"+j+" "+fMes);
                                                tTx="oSet-"+oCount+" Failed Rev="+(rSt+1)+"/"+(rEn+1)+" Probe="+(pSt+1)+"/"+(pEn+1)+" Fwd="+(fSt+1)+"/"+(fEn+1)+" "+fMes;
                                                System.out.println(tTx);
                                                outOut.write(tTx+System.getProperty("line.separator"));
                                            }
                                            test=true;
                                        }
                                    }
                                    //Reverse (3'-5') of the Probe
                                    //for(int k=i+1;k<j-pLens[1];k++) {//old way - with no probe overlaps
                                    for(int k=rSt;k<j-pLens[1];k++) { 
                                        
                                        if(scores[k][4][0]<bestProbe) {
                                            bestProbe=scores[k][4][0];
                                            bestProbeSt=k+pLens[4]+1;
                                            bestProbeEn=k;
                                            bestProbePer=100-(double)scores[k][4][0]/(double)pLens[4]*100;
                                        
                                        }
                                        
                                        pSt=k+pLens[4]+1;
                                        pEn=k;

                                        if(scores[k][4][13]==0) {
                                            if(oPass) {
                                                oCount++;
                                                //System.out.println("S"+oCount+" RT-PCR Success Rev="+(i-pLens[2]+1)+"/"+i+" Probe="+k+"/"+(k-pLens[4]+1)+" Fwd="+(j+pLens[3]-1)+"/"+j);

                                                double cts[]=calculateCT(fEn,3,pEn,4,rEn,2);
                                                tTx="oSet-"+oCount+" RT-PCR Success Rev="+(rSt+1)+"/"+(rEn+1)+" Probe="+(pSt+1)+"/"+(pEn+1)+" Fwd="+(fSt+1)+"/"+(fEn+1);
                                                System.out.println(tTx);
                                                outOut.write(tTx+System.getProperty("line.separator"));
                                                
                                                tTx="deltaCT="+df.format(cts[0])+" [fCT="+df.format(cts[1])+" pCT="+df.format(cts[2])+" rCT="+df.format(cts[3])+"]";
                                                System.out.println(tTx);
                                                outOut.write(tTx+System.getProperty("line.separator"));
                                                
                                                outCT.write(seqNames[s]+"\t"+df.format(cts[0])+"\t"+df.format(cts[4])+"\t"+"Fwd[CT="+df.format(cts[1])+"]="+(fSt+1)+"/"+(fEn+1)+"\tProbe=[CT="+df.format(cts[2])+"]="+(pSt+1)+"/"+(pEn+1)+"\tRev=[CT="+df.format(cts[3])+"]="+(rSt+1)+"/"+(rEn+1)+"\t"+(fSt-rSt)+System.getProperty("line.separator"));
                                                success=true;
                                            }
                                            else {
                                                oCount++;
                                                //System.out.println("S"+oCount+" Failed Rev="+(i-pLens[2]+1)+"/"+i+" Probe="+k+"/"+(k-pLens[1]+1)+" Fwd="+(j+pLens[3]-1)+"/"+j+" "+fMes);
                                                tTx="oSet-"+oCount+" Failed Rev="+(rSt+1)+"/"+(rEn+1)+" Probe="+(pSt+1)+"/"+(pEn+1)+" Fwd="+(fSt+1)+"/"+(fEn+1)+" "+fMes;
                                                System.out.println(tTx);
                                                outOut.write(tTx+System.getProperty("line.separator"));
                                            }
                                            test=true;
                                        }
                                    }

                                    //if no probe site found
                                    if(!test) {
                                       oPass=false;
                                       fMes+="(Rule 3) Failed as without [effective] probe "; 

                                        if(!fMes.equals("")) {
                                            oCount++;
                                            //System.out.println("S"+oCount+" Failed Rev="+(i-pLens[0]+1)+"/"+i+" Fwd="+(j+pLens[2]-1)+"/"+j+" "+fMes);
                                            tTx="oSet-"+oCount+" Failed Rev="+(rSt+1)+"/"+(rEn+1)+" Fwd="+(fSt+1)+"/"+(fEn+1)+" BestFailedProbe="+df.format(bestProbePer)+"% "+(bestProbeSt+1)+"/"+(bestProbeEn+1)+" "+fMes;
                                        
                                            System.out.println(tTx);
                                            outOut.write(tTx+System.getProperty("line.separator"));
                                        }
                                    }
                                }//fwd primer good
                            }//(rev) fwd primer loop
                        }//rev primer good
                    }//(fwd) rev primer loop

                    System.out.println();
                    outOut.write(System.getProperty("line.separator"));

                    if(seqOut) {
                        FileWriter fstream = new FileWriter(outFilename);
                        BufferedWriter out = new BufferedWriter(fstream);
                    
                        if(fullOut)
                            out.write(header);
                        else
                            out.write("Position\t53-Fwd%\t53-Fwd1-4\t53-Fwd\t53-Probe%\t53-Probe1-4\t53-Probe\t53-Rev%\t53-Rev1-4\t53-Rev\t35-Fwd%\t35-Fwd1-4\t35-Fwd\t35-Probe%\t35-Probe1-4\t35-Probe\t35-Rev%\t35-Rev1-4\t35-Rev"+System.getProperty("line.separator"));

                        for(int i=0;i<scores.length;i++) {

                            out.write((i+1)+"\t");//genome position

                            if(!fullOut) {
                                for(int j=0;j<scores[i].length;j++) {
                                    double per=100-(double)scores[i][j][0]/pLens[j]*100;
                                    out.write(df.format(per)+"\t");//% match
                                    out.write(scores[i][j][11]+"\t");//muts 1-4
                                    out.write(scores[i][j][13]+"\t");//successful
                                }
                            }
                            else {
                                for(int j=0;j<scores[i].length;j++) {
                                    for(int k=0;k<scores[i][j].length;k++) {
                                       out.write(scores[i][j][k]+"\t");
                                    }
                                }
                            }

                            out.write(System.getProperty("line.separator"));
                        }
                        
                        out.close();
                    }
                }//matrix outputFile try
                
                catch (Exception e) {
                    e.printStackTrace();
                    System.err.println("Error: " + e.getMessage());
                }
                
                if(!success) {
                    noHits+=seqNames[s]+"\tNOHIT"+System.getProperty("line.separator");
                }
            }//end of comSeq loop s - evaluate each seq

            outOut.close();
            
            outCT.write(noHits);
            outCT.close();
        }//cts outFile try
        catch (Exception e) {
            e.printStackTrace();
            System.err.println("Error: " + e.getMessage());
        }
        System.out.println("GoPrime finished...Exiting");
    }
    
    private static int checkMutation(char pBase, char tBase) {
        
        int mut=0;
        int pSel=0, tSel=0;
       
        for(int i=0;i<compsBase.length;i++) {
            if(pBase==compsBase[i]) {
                pSel=i;
                break;
            }
        }
        
        for(int i=0;i<compsBase.length;i++) {
            if(tBase==compsBase[i]) {
                tSel=i;
                break;
            }
        }
        
        //check if it is a match
        boolean test=false;
        for(int i=0;i<comps[pSel].length;i++) {
            if(comps[pSel][i]==tBase) {
                test=true;
                break;
            }
        }
        
        //if not - check if the mutation is a transition or transversion
        if(!test) {
            //New - we have to revcomp the primer base
            //This is what we expect in the template
            //Compare this to what we observe to determine if it is mutated or not
            pBase=compBase(pBase);
            
            for(int i=0;i<compsBase.length;i++) {
                if(pBase==compsBase[i]) {
                    pSel=i;
                    break;
                }
            }
            
            boolean transi=false;
            for(int i=1;i<muts[pSel].length;i++) {
                if(transi)
                    break;
                
                for(int j=1;j<muts[tSel].length;j++) {
                    if((muts[pSel][i]==ts[0][0] & muts[tSel][j]==ts[0][1])
                      |(muts[pSel][i]==ts[0][1] & muts[tSel][j]==ts[0][0])
                      |(muts[pSel][i]==ts[1][0] & muts[tSel][j]==ts[1][1])
                      |(muts[pSel][i]==ts[1][1] & muts[tSel][j]==ts[1][0])) {
                        mut=1;
                        transi=true;
                        break;
                    }
                }
            }
            if(!transi) {
                mut=2;
            }

        }
        
        return mut;
    }
    
    private static String reverseComplement(String inSeq) {
        
        String rev=reverse(inSeq);
        String comp=complement(rev);

        return comp;
    }
    
    private static String reverse(String inSeq) {
        
        String rev="";
            
        //Reverse
        for(int i=inSeq.length()-1;i>=0;i--) {
            rev=rev+inSeq.charAt(i);
        }

        return rev;
    }
    
    
    private static char compBase(char inBase) {
        
        //Complement
        char comp='N';
        
        if(inBase=='A') {
            comp='T';
        }
        else if(inBase=='C') {
            comp='G';
        }
        else if(inBase=='G') {
            comp='C';
        }
        else if(inBase=='T') {
            comp='A';
        }
        else {
            boolean test=false;
            for(int j=0;j<revCodes.length;j++) {
                if(revCodes[j][0]==inBase) {
                    comp=revCodes[j][1];
                    test=true;
                    break;
                }
            }
            if(!test) {
                System.out.println("ERROR - unrecognised character during complementation "+inBase);
            }
        }
        
        return comp;
    }
        
    
    private static String complement(String inSeq) {
        
        //Complement
        String comp="";
        for(int i=0;i<inSeq.length();i++) {
            if(inSeq.charAt(i)=='A') {
                comp=comp+"T";
            }
            else if(inSeq.charAt(i)=='C') {
                comp=comp+"G";
            }
            else if(inSeq.charAt(i)=='G') {
                comp=comp+"C";
            }
            else if(inSeq.charAt(i)=='T') {
                comp=comp+"A";
            }
            else {
                boolean test=false;
                for(int j=0;j<revCodes.length;j++) {
                    if(revCodes[j][0]==inSeq.charAt(i)) {
                        comp=comp+revCodes[j][1];
                        test=true;
                        break;
                    }
                }
                if(!test) {
                    System.out.println("ERROR - unrecognised character during complementation "+inSeq.charAt(i));
                }
            }
        }

        return comp;
    }

    
    private static String checkSeqs(String name, String inSeq) {
        
        String outSeq="";
        String bad="";
        int b=0;
        
        //Switch Us to Ts, Remove GAPs
        outSeq=inSeq.toUpperCase().replace("U", "T").replace("-", "");
        
        for(int i=0;i<outSeq.length();i++) {
            boolean test=false;
            
            for(int j=0;j<allow.length;j++) {
                if(outSeq.charAt(i)==allow[j]) {
                    test=true;
                    break;
                }
            }
            if(!test) {
                bad+=outSeq.charAt(i)+" ";
                b++;
            }
        }
        
        if(b>0) {
            System.out.println(name+" had "+b+" bad seq bases = "+bad+" original seq = "+inSeq);
            System.out.println("Exiting...");
            System.exit(0);
        }

        return outSeq;
    }
    
    private static void setUpBases() {
        
        head="TotMis\tTotTs\tTotTv\tNt1Ts\tNt1Tv\tNt2Ts\tNt2Tv\tNt3Ts\tNt3Tv\tNt4Ts\tNt4Tv\t1-4NtMis\tMisPass\tPass\tNs\t5-1-4NtMis";
        header="";
        header="Pos\t"+"[Fwd]"+head+"\t[Probe]"+head+"\t[Rev]"+head+"\t[Fwd]"+head+"\t[Probe]"+head+"\t[Rev]"+head+System.getProperty("line.separator");
        
        //Allowed DNA characters
        allow=new char[16];
        allow[0]='A';
        allow[1]='C';
        allow[2]='G';
        allow[3]='T';//1st thing we do is flip Us to Ts - so Us not allowed technically - so not in list
        allow[4]='M';
        allow[5]='R';
        allow[6]='W';
        allow[7]='S';
        allow[8]='Y';
        allow[9]='K';
        allow[10]='V';
        allow[11]='H';
        allow[12]='D';
        allow[13]='B';
        allow[14]='N';
        allow[15]='-';//left in for future alignment stuff - technically remove gaps at start currently
        
        //Bases and revComps
        bases=new char[4][2];
        bases[0][0]='A';bases[0][1]='T';
        bases[1][0]='C';bases[1][1]='G';
        bases[2][0]='G';bases[2][1]='C';
        bases[3][0]='T';bases[3][1]='A';
        
        //The Reverse Complement of IUPAC Ambiguity codes
        revCodes=new char[12][2];
        revCodes[0][0]='M';revCodes[0][1]='K';
        revCodes[1][0]='R';revCodes[1][1]='Y';
        revCodes[2][0]='W';revCodes[2][1]='W';
        revCodes[3][0]='S';revCodes[3][1]='S';
        revCodes[4][0]='Y';revCodes[4][1]='R';
        revCodes[5][0]='K';revCodes[5][1]='M';
        revCodes[6][0]='V';revCodes[6][1]='B';
        revCodes[7][0]='H';revCodes[7][1]='D';
        revCodes[8][0]='D';revCodes[8][1]='H';
        revCodes[9][0]='B';revCodes[9][1]='V';
        revCodes[10][0]='N';revCodes[10][1]='N';
        revCodes[11][0]='-';revCodes[11][1]='-';
        
        //This is what IUPAC amigiuity codes correspond to
        amb=new char[11][];
        amb[0]=new char[3];
        amb[0][0]='M';amb[0][1]='A';amb[0][2]='C';
        amb[1]=new char[3];
        amb[1][0]='R';amb[1][1]='A';amb[1][2]='G';
        amb[2]=new char[3];
        amb[2][0]='W';amb[2][1]='A';amb[2][2]='T';
        amb[3]=new char[3];
        amb[3][0]='S';amb[3][1]='C';amb[3][2]='G';
        amb[4]=new char[3];
        amb[4][0]='Y';amb[4][1]='C';amb[4][2]='T';
        amb[5]=new char[3];
        amb[5][0]='K';amb[5][1]='G';amb[5][2]='T';
        amb[6]=new char[4];
        amb[6][0]='V';amb[6][1]='A';amb[6][2]='C';amb[6][3]='G';
        amb[7]=new char[4];
        amb[7][0]='H';amb[7][1]='A';amb[7][2]='C';amb[7][3]='T';
        amb[8]=new char[4];
        amb[8][0]='D';amb[8][1]='A';amb[8][2]='G';amb[8][3]='T';
        amb[9]=new char[4];
        amb[9][0]='B';amb[9][1]='C';amb[9][2]='G';amb[9][3]='T';
        amb[10]=new char[5];
        amb[10][0]='N';amb[10][1]='A';amb[10][2]='C';amb[10][3]='G';amb[10][4]='T';
        
        int cCount=0;
        String cText="";
        char bComps[][]=new char[4][];
        
        //loop to get all the reverse complements of each base including ambiguity codes
        for(int b=0;b<bases.length;b++) {
            cCount=1;
            cText=""+bases[b][1];
            
            for(int i=0;i<amb.length;i++) {
                for(int j=1;j<amb[i].length;j++) {
                    if(amb[i][j]==bases[b][1]) {//NB - [b][1] - searching for the reverse complement of the base [b][0]
                        cText+=amb[i][0];
                        cCount++;
                    }
                }
            }
            
            //for each base b - store all possible revcomps bases/amiguities
            bComps[b]=new char[cText.length()];
            for(int i=0;i<cText.length();i++) {
                bComps[b][i]=cText.charAt(i);
            }
        }
       
        char aComps[][]=new char[amb.length][];
        
        //loop to find all the reverse complements of each amiguity code
        for(int a=0;a<amb.length;a++) {
            cCount=0;
            cText="";
            
            //loop through all the bases the amiguity stands for
            for(int j=1;j<amb[a].length;j++) {
                //for each base that the amiguity represents - find all its complements (incl ambiguity codes)
                for(int k=0;k<bases.length;k++) {
                    if(amb[a][j]==bases[k][0]) {
                        //bComps[k] stores all the rev comps of bases[k]
                        for(int l=0;l<bComps[k].length;l++) {
                            //check the base/code is not used already in cText
                            boolean test=false;
                            for(int q=0;q<cText.length();q++) {
                                if(bComps[k][l]==cText.charAt(q)) {
                                    test=true;
                                }
                            }
                            //if not add it
                            if(!test) {
                               cText+=bComps[k][l]; 
                               cCount++;
                            }
                        }
                        break;
                    }
                }//end of bases[k] loop
            }//end of amb[i][j] loop
            
            aComps[a]=new char[cText.length()];
            for(int j=0;j<cText.length();j++) {
                aComps[a][j]=cText.charAt(j);
            }
        }//end of amb[a] loop
        
        //Basically compsBase stores the actual base or ambiguity code
        //Whilst comps stores all the rev comps [bases and ambiguity] for the compsBase
        comps=new char[bComps.length+aComps.length][];
        for(int i=0;i<bComps.length;i++) {
            comps[i]=bComps[i];
        }
        for(int i=0;i<aComps.length;i++) {
            comps[i+bComps.length]=aComps[i];
        }
        
        compsBase=new char[comps.length];
        for(int i=0;i<bases.length;i++) {
            compsBase[i]=bases[i][0];
        }
        for(int i=0;i<amb.length;i++) {
            compsBase[i+bases.length]=amb[i][0];
        }
        
        
        //Transitions - no need to store transversion - if it is not a Ts it is a Tv
        ts=new char[2][2];
        ts[0][0]='A';ts[0][1]='G';//AG -> purines
        ts[1][0]='C';ts[1][1]='T';//CT -> pyrimidines
        
        //this is essentially a merger of bases and amb
        //single array to store what each base/amiguity represents
        //the bases represent themselves!
        //used in check mutations for checking TsTv
        muts=new char[bases.length+amb.length][];
        muts[0]=new char[2];muts[1]=new char[2];muts[2]=new char[2];muts[3]=new char[2];
        muts[0][0]='A';muts[0][1]='A';
        muts[1][0]='C';muts[1][1]='C';
        muts[2][0]='G';muts[2][1]='G';
        muts[3][0]='T';muts[3][1]='T';
        for(int i=0;i<amb.length;i++) {
            muts[i+bases.length]=amb[i];//amb[i] is an array
        }
        
    }//end of setUpBases()
     
    public static double[] calculateCT(int fPos, int fDir, int pPos, int pDir, int rPos, int rDir) {
        
        /*
        f = forward, p = probe, r = reverse
        fPos pPos rPos store the sequence position the 3' end of the primer binds - where the data is stored
        fDir pDir rDir store the direction [0-6] - where the data is in 2nd array dimension
        
        scores[seqLenth][5'-3'(FwdProbRev)3'-5'(FwdProbeRev)][14]
        0 total mismatchs
        1 total transition
        2 total transversion
        3 nt 1 transition mismatch
        4 nt 1 transversion mismatch
        5 nt 2 transition mismatch
        6 nt 2 transversion mismatch
        7 nt 3 transition mismatch
        8 nt 3 transversion mismatch
        9 nt 4 transition mismatch
        10 nt 4 transversion mismatch
        11 mismatches in 1-4
        12 >=82.05% (primers) or >=85% (probe) matches
        13 candidate
        14 Ns
        15 mistaches at 5' 1-4 (for probe)
        */ 
        
        //should replace all below with an intelligent loop
        
        double ct=0, fCt=0, rCt=0, pCt=0;
        int fCount=0, pCount=0, rCount=0, pCount5=0;
        int f2=0, r2=0, b2=0;
        //forward
        if(scores[fPos][fDir][3]>0) {
            fCt+=primerCts[1]*scores[fPos][fDir][3];//nt 1 Ts
            fCount++;
        }
        if(scores[fPos][fDir][4]>0) {
            fCt+=primerCts[2]*scores[fPos][fDir][4];//nt 1 Tv
            fCount++;
        }
        if(scores[fPos][fDir][5]>0) {
            fCt+=primerCts[3]*scores[fPos][fDir][5];//nt 2 Ts
            fCount++;
        }
        if(scores[fPos][fDir][6]>0) {
            fCt+=primerCts[4]*scores[fPos][fDir][6];//nt 2 Tv
            fCount++;
        }
        if(scores[fPos][fDir][7]>0) {
            if(scores[fPos][fDir][9]>0) {
                fCt+=primerCts2[5];
                f2++;
            }
            else
                fCt+=primerCts[5]*scores[fPos][fDir][7];//nt 3 Ts
            
            fCount++;
        }
        if(scores[fPos][fDir][8]>0) {
            if(scores[fPos][fDir][10]>0) {
                fCt+=primerCts2[6];
                f2++;
            }
            else
                fCt+=primerCts[6]*scores[fPos][fDir][8];//nt 3 Tv
            
            fCount++;
        }
        if(scores[fPos][fDir][9]>0) {
            if(scores[fPos][fDir][7]==0)
                fCt+=primerCts[7]*scores[fPos][fDir][9];//nt 4 Ts
            
            fCount++;
        }
        if(scores[fPos][fDir][10]>0) {
            if(scores[fPos][fDir][8]==0)
                fCt+=primerCts[8]*scores[fPos][fDir][10];//nt 4 Tv
            
            fCount++;
        }
        if(scores[fPos][fDir][0]>0) {
            fCt+=(double)scores[fPos][fDir][0]/(double)mLen*100*(double)primerCts[0];//%
            //fCount++;
        }
        
        //reverse
        if(scores[rPos][rDir][3]>0) {
            if(scores[fPos][fDir][3]>0) {
                rCt+=primerCts2[1]/2;//nt 1 Ts
                fCt+=(primerCts2[1]/2-primerCts[1]);
                b2++;
            }
            else
                rCt+=primerCts[1]*scores[rPos][rDir][3];//nt 1 Ts
            
            rCount++;
        }
        if(scores[rPos][rDir][4]>0) {
            if(scores[fPos][fDir][4]>0) {
                rCt+=primerCts2[2]/2;//nt 1 Tv
                fCt+=(primerCts2[2]/2-primerCts[2]);//nt 1 Tv
                b2++;
            }
            else
                rCt+=primerCts[2]*scores[rPos][rDir][4];//nt 1 Tv
            
            rCount++;
        }
        if(scores[rPos][rDir][5]>0) {
            if(scores[fPos][fDir][5]>0) {
                rCt+=primerCts2[3]/2;//nt 2 Ts
                fCt+=(primerCts2[3]/2-primerCts[3]);//nt 2 Ts
                b2++;
            }
            else
                rCt+=primerCts[3]*scores[rPos][rDir][5];//nt 2 Ts
            
            rCount++;
        }
        if(scores[rPos][rDir][6]>0) {
            if(scores[fPos][fDir][6]>0) {
                rCt+=primerCts2[4]/2;//nt 2 Tv
                fCt+=(primerCts2[4]/2-primerCts[4]);//nt 2 Tv
                b2++;
            }
            else
                rCt+=primerCts[4]*scores[rPos][rDir][6];//nt 2 Tv
            
            rCount++;
        }
        
        if(scores[rPos][rDir][7]>0) {
            if(scores[fPos][fDir][7]>0 | scores[fPos][fDir][9]>0) {//if fwd primer also has
                rCt+=primerCts2[5]/2;//nt 3-4 Ts
                fCt+=(primerCts2[5]/2-primerCts[5]);//nt 3-4 Ts
                b2++;
            }
            //As we have max 2 mutations in 1-4 between primers - this if/else works
            else if(scores[rPos][rDir][9]>0) {//if pos 4 in rev also has
                rCt+=primerCts2[5];
                r2++;
            }
            else
                rCt+=primerCts[5]*scores[rPos][rDir][7];//nt 3 Ts
            
            rCount++;
        }
        if(scores[rPos][rDir][8]>0) {
            if(scores[fPos][fDir][8]>0 | scores[fPos][fDir][10]>0) {
                rCt+=primerCts2[6]/2;//nt 3-4 Tv
                fCt+=(primerCts2[6]/2-primerCts[6]);//nt 3-4 Tv
                b2++;
            }
            else if(scores[rPos][rDir][10]>0){
                rCt+=primerCts2[6];
                r2++;
            }
            else
                rCt+=primerCts[6]*scores[rPos][rDir][8];//nt 3 Tv
            
            rCount++;
        }
        if(scores[rPos][rDir][9]>0) {
            if(scores[fPos][fDir][7]>0 | scores[fPos][fDir][9]>0) {
                rCt+=primerCts2[7]/2;//nt 3-4 Ts
                fCt+=primerCts2[7]/2-primerCts[7];//nt 3-4 Ts
                b2++;
            }
            else if(scores[rPos][rDir][7]>0) {
                //do nothing as already added the combined affect above
            }
            else 
                rCt+=primerCts[7]*scores[rPos][rDir][9];//nt 4 Ts
            
            rCount++;
        }
        if(scores[rPos][rDir][10]>0) {
            if(scores[fPos][fDir][8]>0 | scores[fPos][fDir][10]>0) {
                rCt+=primerCts2[8]/2;//nt 3-4 Tv
                fCt+=primerCts2[8]/2-primerCts[8];//nt 3-4 Tv
                b2++;
            }
            else if(scores[rPos][rDir][8]>0) {
                //do nothing as already added the combined affect above
            }
            else
                rCt+=primerCts[8]*scores[rPos][rDir][10];//nt 4 Tv
            
            rCount++;
        }
        if(scores[rPos][rDir][0]>0) {
            rCt+=(double)scores[rPos][rDir][0]/(double)mLen*100*(double)primerCts[0];//%
            //rCount++;
        }
        
        //Prob want to report/record 3' and 5' differently
        //probe
        if(scores[pPos][pDir][11]>0) {//total mismatches in 3' 1-4
            pCt+=probeCts[1]*scores[pPos][pDir][11];
            pCount+=scores[pPos][pDir][11];//
        }
        if(scores[pPos][pDir][15]>0) {//total mismatches in 5' 1-4
            pCt+=probeCts[5]*scores[pPos][pDir][15];
            pCount5+=scores[pPos][pDir][15];
        }
        
        /*
        //Old way
        if(scores[pPos][pDir][3]>0) {
            pCt+=probeCts[1]*scores[pPos][pDir][3];//nt 1 Ts
            pCount++;
        }

        if(scores[pPos][pDir][4]>0) {
            pCt+=probeCts[2]*scores[pPos][pDir][4];//nt 1 Tv
            pCount++;
        }

        int ts24=scores[pPos][pDir][5]+scores[pPos][pDir][7]+scores[pPos][pDir][9];
        int tv24=scores[pPos][pDir][6]+scores[pPos][pDir][8]+scores[pPos][pDir][10];
        
        if(ts24>0) {
            pCt+=probeCts[3]*ts24;
            pCount++;
        }

        if(tv24>0) {
            pCt+=probeCts[4]*tv24;
            pCount++;
        }
        */

        if(scores[pPos][pDir][0]>0) {
           pCt+=(double)scores[pPos][pDir][0]/(double)pLens[pDir]*100*(double)probeCts[0];//%
           //pCount++;
        }

        if(fCount>1 & f2==0)
            fCt+=intraPrimer;
        if(rCount>1 & r2==0)
            rCt+=intraPrimer;
        if(pCount>1)
            pCt+=intraProbe;

        if(pCt<0)
            pCt=0;
        
        ct=fCt+pCt+rCt;
        
        if(fCount>0 & rCount>0 & b2==0)
            ct+=interPrimer;
        if(pCount>0 & (fCount>0 | rCount>0))
            ct+=interProbe;
                
        double cts[]=new double[5];
        cts[0]=ct;
        cts[1]=fCt;
        cts[2]=pCt;
        cts[3]=rCt;
        cts[4]=ct/deltaLod;
        
        //System.out.println("CT="+df.format(cts[0])+" [fCT="+df.format(cts[1])+" pCT="+df.format(cts[2])+" rCT="+df.format(cts[3])+"]");
        
        return cts;
    }
    
    public static void setUpCts() {
        
        minPm=1-0.8205;//maximum Primer mismatch % based on combined F and R length
        minPmI=1-0.7;//maximum INDIVIDUAL Primer mismacth % for an individual primer
        minPb=1-0.85;//maximum Probe mismtach %
         
        primerCts=new double[9];//this could do down to 7 now
        primerCts2=new double[9];//and this
         
        primerCts[0]=0.86313;//per % - was 1.05177
        primerCts[1]=1.63934;//nt 1 Ts - was 3.351592
        primerCts[2]=4.0097;//nt 1 Tv - was 5.762987
        primerCts[3]=1.02577;//nt 2 Ts - was 2.460959
        primerCts[4]=3.64175;//nt 2 Tv - was 3.719203
        primerCts[5]=1.18575;//nt 3-4 Ts - was 0.981312 for just nt 3
        primerCts[6]=2.91193;//nt 3-4 Tv - was 2.007043 for just nt 3
        primerCts[7]=1.18575;//nt 4 Ts - was 1.203353 for just nt 4
        primerCts[8]=2.91193;//nt 4 Tv - was 2.42565 for just nt 4
        
        primerCts2[0]=0;//% doesnt have a differnet value for combined
        primerCts2[1]=5.06991;//2 Ts at nt 1
        primerCts2[2]=8.89553;//2 Tv at nt 1
        primerCts2[3]=3.51227;//2 Ts at nt 2
        primerCts2[4]=6.31407;//2 Tv at nt 2
        primerCts2[5]=1.7551;//2 Ts at nt 3-4
        primerCts2[6]=4.85546;//2 Tv at nt 3-4
        primerCts2[7]=1.7551;//same as [5] - just including for completeness
        primerCts2[8]=4.85546;//same as [6] - just including for completeness

        probeCts=new double[6];
        probeCts[0]=0.35490;//per % - was 0.282996
        probeCts[1]=0.40779;//was 1.767975;//nt 1 Ts - 1-4 TsTv now
        probeCts[2]=0.40779;//was 2.158492;//nt 1 Tv - 1-4 TsTv now
        probeCts[3]=0.40779;//was 1.6015;//nt 2-4 Ts - 1-4 TsTv now
        probeCts[4]=0.40779;//was 2.125075;//nt 2-4 Tv - 1-4 TsTv now
        probeCts[5]=-1.04937;//New 5' 1-4 nt
       
        intraPrimer=0.27264;//was -0.9923;
        interPrimer=0.27264;//was -0.8052;
        intraProbe=-0.83874;//not sure if this still exists
        interProbe=0.42976;//was 1.6729;//inter as in probe with primer
        
        deltaLod=3.43527;
        
        max14nt=2;//maximum mutatins in 1-4 nts at 3' end of primer
        
        fullOut=true;
        seqOut=false;
         
    }//end of setUpCts
    
}
