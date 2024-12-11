
import com.github.sh0nk.matplotlib4j.Plot;
import com.github.sh0nk.matplotlib4j.PythonExecutionException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

public class Test
{
    public static void HashAdd(HashMap<Double,Double> h, Double key, Double val)
    {
        if (h.containsKey(key))
        {
            h.put(key, h.get(key) + val);
        }
        else
        {
            h.put(key, val);
        }
    }
    public static void Run(int testsToDo, boolean decay, double mutationAdjustement, int queryLength, ProbabilisticDatabase d, double gapProbabilitySequence, int iterations, int wordSize, double indexingThreshold, double ungappedExtensionThreshold,double delta, int gapPenalty, double safetyZero)
    {
        Random rand = new Random();
        String[] nucleotide = {"A","G","C","T"};
        ArrayList<Integer> iterList = new ArrayList<>();
        ArrayList<Integer> indexMatchList = new ArrayList<>();
        ArrayList<Integer> ungappedMatchList = new ArrayList<>();
        ArrayList<Integer> positiveList = new ArrayList<>();
        String nucleotideSequence = "";
        //Determine start and stop of sequence
        int Start = 20000;//(int)(Math.random() * (d.size() - queryLength));
        int Stop = Start + queryLength;
        System.out.println(Start + " " + Stop);
        HashMap<Integer, Integer> index = new HashMap<Integer, Integer>();
        //Create sequence assuming a probability of 1.
        for (int n = Start;n<Stop;n++) {
            //Add deletion / insertion
            if (Math.random() < gapProbabilitySequence) {
                //Insertion
                if (Math.random() < 0.5) {
                    nucleotideSequence += d.getNucleotide(n);
                    index.put(nucleotideSequence.length() - 1, n);
                    nucleotideSequence += nucleotide[rand.nextInt(4)];
                }
                //Else deletion, append nothing
                else {
                    nucleotideSequence += d.getNucleotide(n);
                    index.put(nucleotideSequence.length() - 1, n);
                }
            }
            else
            {
                nucleotideSequence += d.getNucleotide(n);
                index.put(nucleotideSequence.length() - 1, n);
            }
        }
        Sequence s = new Sequence(nucleotideSequence);
        //Run tests
        int positive=0;
        int falsePositive=0;
        Index i = new Index(d, wordSize);
        for(int testNumber=0;testNumber<testsToDo;testNumber++)
        {
            positive=0;
            falsePositive=0;
            System.out.println("Running Test " + testNumber + " ***************************************************************************************");
            String[] out = BLAST.ProbabiliticBlast(s, i, d, wordSize, indexingThreshold, ungappedExtensionThreshold, delta, iterations, gapPenalty, safetyZero);
            if (Double.valueOf(out[0]) >= Start && Double.valueOf(out[1]) <= Stop)
            {
                positive += 1;
            }
            else
            {
                falsePositive +=1;
            }
            System.out.println("Word matches: " + out[5] + " Ungappped matches: " + out[6]);
            System.out.println("P: " + positive + " F: " + falsePositive);
            if (!out.equals(new String[]{"0", "0", "0", "0", "0", "0", "0"}))
            {
                System.out.println("Final Sequence Match:");
                System.out.println(out[2]);
                System.out.println(out[3]);
                positiveList.add(positive);
                indexMatchList.add(Integer.valueOf(out[5]));
                ungappedMatchList.add(Integer.valueOf(out[6]));
                iterList.add(testNumber);
            }
            //Degrade sequence
            if (decay) {
                if (testNumber + 1 < iterations) {
                    String nucleotideSequence2 = "";
                    for (int n = 0; n < nucleotideSequence.length(); n++) {
                        if (index.containsKey(n)) {
                            nucleotideSequence2 += (Math.random() < (d.getProbability(index.get(n))) - mutationAdjustement) ? String.valueOf(nucleotideSequence.charAt(n)) : nucleotide[rand.nextInt(4)];
                        } else {
                            nucleotideSequence2 += String.valueOf(nucleotideSequence.charAt(n));
                        }
                    }
                    nucleotideSequence = nucleotideSequence2;
                    s = new Sequence(nucleotideSequence);
                }
            }
        }
        Plot plt1 = Plot.create();
        Plot plt2 = Plot.create();
        Plot plt3 = Plot.create();
        plt1.plot().add(iterList, indexMatchList);
        plt2.plot().add(iterList,ungappedMatchList);
        plt3.plot().add(iterList,positiveList);
        plt1.title("Index matches");
        plt2.title("Ungapped matches");
        plt3.title("Final Sequence Correctness");
        plt1.xlabel("Mutation Rounds");
        plt2.xlabel("Mutation Rounds");
        plt3.xlabel("Mutation Rounds");
        plt1.ylabel("Found Sequences");
        plt2.ylabel("Found Sequences");
        plt3.ylabel("Correct (1), Incorrect (0)");
        //plt1.xticks(iterList);
        //plt2.xticks(iterList);
        try {
            plt1.show();
            plt2.show();
            plt3.show();
        }catch (PythonExecutionException e)
        {
            System.out.println("Failed to make plot \uD83D\uDE31");
        }
        catch (java.io.IOException e)
        {
            System.out.println("Failed to make plot \uD83D\uDE31");
        }
        catch (NullPointerException e)
        {
            System.out.println("No match \uD83D\uDE31");
        }
    }
    public static void ReadAndRun(File f, ProbabilisticDatabase d, int iterations, int wordSize, double indexingThreshold, double ungappedExtensionThreshold, double delta, int gapPenalty, double safetyZero)
    {
        //Create list of sequences from file
        ArrayList<Sequence> sequencesToTest = new ArrayList<>();
        ArrayList<Integer> sequencesToTestIndex = new ArrayList<>();
        ArrayList<Double> sequencesToTestProb = new ArrayList<>();
        try {
            Scanner testReader = new Scanner(f);
            while (testReader.hasNextLine())
            {
                sequencesToTestProb.add(Double.parseDouble(testReader.nextLine()));
                sequencesToTestIndex.add(Integer.valueOf(testReader.nextLine()));
                sequencesToTest.add(new Sequence(testReader.nextLine()));
            }
            testReader.close();
        }catch (FileNotFoundException e)
        {
            System.out.println("Could not find file \uD83D\uDE31");
        }
        //Run test for every sequence
        HashMap<Double, Double> AccuracyPerProb = new HashMap<>();
        HashMap<Double, Double> NumberOfSeqPerProb = new HashMap<>();
        HashMap<Double,Double> IndexPerProb = new HashMap<>();
        HashMap<Double,Double> UngappedPerProb = new HashMap<>();
        Index i = new Index(d, wordSize);
        for(int testNumber=0; testNumber<sequencesToTest.size()/3;testNumber++)
        {
            Sequence s = sequencesToTest.get(testNumber);
            int Start = sequencesToTestIndex.get(testNumber);
            int Stop = Start + s.size();
            System.out.println("Test " + testNumber + " Start: " + Start + " Stop: " + Stop + " Prob: " + sequencesToTestProb.get(testNumber));
            //Run tests
            double positive=0;
            double falsePositive=0;
            String[] out = BLAST.ProbabiliticBlast(s, i, d, wordSize, indexingThreshold, ungappedExtensionThreshold, delta, iterations, gapPenalty, safetyZero);
            if ((Double.valueOf(out[0]) >= (Double.valueOf(out[0]) - s.size() * sequencesToTestProb.get(testNumber))) && (Double.valueOf(out[0]) <= (Double.valueOf(out[0]) + s.size() * sequencesToTestProb.get(testNumber))) && !Arrays.equals(out,new String[]{"0", "0", "0", "0", "0", "0", "0"}))
            {
                positive += 1;
            }
            else
            {
                falsePositive +=1;
            }
            System.out.println("Word matches: " + out[5] + " Ungappped matches: " + out[6]);
            System.out.println("P: " + positive + " F: " + falsePositive);
            if (!Arrays.equals(out,new String[]{"0", "0", "0", "0", "0", "0", "0"})) {
                System.out.println("Final Sequence Match:");
                System.out.println(out[2]);
                System.out.println(out[3]);
                System.out.println(out[4]);
            }
            HashAdd(AccuracyPerProb, sequencesToTestProb.get(testNumber),positive);
            HashAdd(NumberOfSeqPerProb, sequencesToTestProb.get(testNumber),1.);
            HashAdd(IndexPerProb, sequencesToTestProb.get(testNumber),Double.valueOf(out[5]));
            HashAdd(UngappedPerProb, sequencesToTestProb.get(testNumber),Double.valueOf(out[6]));
        }
        //Calculate average for each
        // Print keys and values
        ArrayList<Double> probValues = new ArrayList<>(NumberOfSeqPerProb.keySet());
        ArrayList<Double> indexValues = new ArrayList<>();
        ArrayList<Double> ungappedValues = new ArrayList<>();
        ArrayList<Double> accuracyValues = new ArrayList<>();
        for (double p : probValues)
        {
            accuracyValues.add(AccuracyPerProb.get(p) / NumberOfSeqPerProb.get(p));
            indexValues.add(IndexPerProb.get(p) / NumberOfSeqPerProb.get(p));
            ungappedValues.add(UngappedPerProb.get(p) / NumberOfSeqPerProb.get(p));
        }

        Plot plt1 = Plot.create();
        Plot plt2 = Plot.create();
        Plot plt3 = Plot.create();
        plt1.plot().add(probValues, accuracyValues);
        plt2.plot().add(probValues,indexValues);
        plt3.plot().add(probValues,ungappedValues);
        plt1.title("Accuracy");
        plt2.title("Index Matches");
        plt3.title("Ungapped Matches");
        plt1.xlabel("Mutation Probability");
        plt2.xlabel("Mutation Probability");
        plt3.xlabel("Mutation Probability");
        plt1.ylabel("Percentage Correct Match");
        plt2.ylabel("Matches");
        plt3.ylabel("Matches");
        try {
            plt1.show();
            plt2.show();
            plt3.show();
        }catch (PythonExecutionException e)
        {
            System.out.println("Failed to make plot \uD83D\uDE31");
        }
        catch (java.io.IOException e)
        {
            System.out.println("Failed to make plot \uD83D\uDE31");
        }
        catch (NullPointerException e)
        {
            System.out.println("No match \uD83D\uDE31");
        }
    }
}
