
import com.github.sh0nk.matplotlib4j.Plot;
import com.github.sh0nk.matplotlib4j.PythonExecutionException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

public class Test
{
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
}
