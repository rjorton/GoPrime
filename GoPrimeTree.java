package goprimetree;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class GoPrimeTree {
    
    public static String cols[];
    public static double colsN[];
    public static int count;
    
    public static String ctsFilename, treeFilename, outFilename;

    public static void main(String[] args) {

        System.out.println("GoPrimeTree started...");
        
        if(args.length==2) {
            ctsFilename=args[0];
            treeFilename=args[1];
        }
        else {
            System.out.println("Incorrect usage - correct usage is java -jar GoPrimeTree.jar CTsFileName TreeFileName ");
            System.exit(1);
        }

        outFilename=treeFilename.substring(0, treeFilename.indexOf(".nex"))+"_col.nex";
        
        System.out.println("CTs filename = "+ctsFilename);
        System.out.println("Tree filename = "+treeFilename);
        System.out.println("Output Tree filename = "+outFilename);
        
        cols=new String[6];
        /*
        cols[0]="[&!color=#bd0026]";
        cols[1]="[&!color=#f03b20]";
        cols[2]="[&!color=#fd8d3c]";
        cols[3]="[&!color=#feb24c]";
        cols[4]="[&!color=#fed976]";
        cols[5]="[&!color=#ffffb2]";
        */

        /*
        cols[0]="[&!color=#1a9850]";
        cols[1]="[&!color=#91cf60]";
        cols[2]="[&!color=#d9ef8b]";
        cols[3]="[&!color=#fee08b]";
        cols[4]="[&!color=#fc8d59]";
        cols[5]="[&!color=#d73027]";
        */
        
        //These are the colours used from ColorBrewer
        cols[0]="[&!color=#1a9850]";
        cols[1]="[&!color=#66bd63]";
        cols[2]="[&!color=#a6d96a]";
        cols[3]="[&!color=#fdae61]";
        cols[4]="[&!color=#f46d43]";
        cols[5]="[&!color=#d73027]";
        
        //These are the increments for CTs
        colsN=new double[cols.length];
        colsN[0]=0;
        colsN[1]=5;
        colsN[2]=10;
        colsN[3]=15;
        colsN[4]=20;
        colsN[5]=25;
        
        count=0;
        
        //Read in the cts File
        File ctsFile = new File(ctsFilename);
        count=0;
        try {
            BufferedReader input =  new BufferedReader(new FileReader(ctsFile));

            try {
                String line = null;

                while (( line = input.readLine()) != null) {
                    if(line.equals(""))
                        continue;
                    
                    String splits[]=line.split("\t");
                    
                    if(!splits[1].equals("NOHIT"))
                        count++;
                }
            }
            finally {
                input.close();
            }
        }
        catch (IOException ex) {
            ex.printStackTrace();
        }
        
        System.out.println(count+" CTs in file "+ctsFile);
        
        //Add in NO-HIT filter
        
        double cts[]=new double[count];
        String seqs[]=new String[count];
        String ctCols[]=new String[count];
        
        count=0;
        String prevSeq="";
        double prevCT=1000000;
        try {
            BufferedReader input =  new BufferedReader(new FileReader(ctsFile));

            try {
                String line = null;

                while (( line = input.readLine()) != null) {
                    if(line.equals(""))
                        continue;
                    
                    String splits[]=line.split("\t");
                    
                    if(!splits[1].equals("NOHIT")) {
                        seqs[count]=splits[0].substring(1, splits[0].length());
                        //seqs[count]=splits[0].substring(1, splits[0].indexOf("."));
                        cts[count]=Double.parseDouble(splits[1]);
                        ctCols[count]=cols[cols.length-1];

                        //Multiple CTs being outputted
                        if(seqs[count].equals(prevSeq)) {
                            System.out.println("Multiple CT found for "+seqs[count]);
                            if(prevCT<cts[count]) {
                                cts[count]=prevCT;
                            }
                            else if(count>0) {
                                cts[count-1]=cts[count];
                            }
                        }

                        prevSeq=seqs[count];
                        prevCT=cts[count];

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
        
        for(int i=0;i<cts.length;i++) {
            for(int j=0;j<colsN.length;j++) {
                if(cts[i]<=colsN[j]) {
                    ctCols[i]=cols[j];
                    break;
                }
            }
        }
        
        
        count=0;
        
        try {
            FileWriter fstream = new FileWriter(outFilename);
            BufferedWriter out = new BufferedWriter(fstream);
        
            File treeFile = new File(treeFilename);
            count=0;
            boolean test=false;
            try {
                BufferedReader input =  new BufferedReader(new FileReader(treeFile));

                try {
                    String line = null;

                    while (( line = input.readLine()) != null) {

                        if(line.indexOf("taxlabels")>=0) {
                            test=true;
                        }
                        else if(test) {
                            if(line.indexOf(";")>=0) {
                                test=false;
                            }
                            else {
                                String seq="";
                                if(line.indexOf("'")>=0) {   
                                    seq=line.substring(line.indexOf("'")+1, line.lastIndexOf("'"));
                                }
                                else {
                                    //if sequence is just a number like >1 >2 - then FigTree does not add ' around it
                                    seq=line.substring(1, line.length());
                                }
                                //System.out.println(seq);
                                for(int i=0;i<seqs.length;i++) {
                                    if(seqs[i].equals(seq)) {
                                        line=line+ctCols[i];
                                        count++;
                                        break;
                                    }

                                }
                            }
                        }

                        out.write(line+System.getProperty("line.separator"));
                    }
                }
                finally {
                    input.close();
                }
            }
            catch (IOException ex) {
                ex.printStackTrace();
            }
               
            out.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.err.println("Error: " + e.getMessage());
        }
        
        System.out.println(count+" colours added to tax in nex file "+treeFilename);
        System.out.println("New tree file = "+outFilename);
        
        System.out.println("GoPrimeTree...exiting");
    }
    
}
