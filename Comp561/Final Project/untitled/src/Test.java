
import com.github.sh0nk.matplotlib4j.Plot;
import com.github.sh0nk.matplotlib4j.PythonExecutionException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

public class Test
{
    public static void Run(int testsToDo, double mutationAdjustement, int queryLength, ProbabilisticDatabase d, double gapProbabilitySequence, int iterations, int wordSize, double indexingThreshold, double ungappedExtensionThreshold,double delta, int gapPenalty)
    {
        Random rand = new Random();
        String[] nucleotide = {"A","G","C","T"};
        ArrayList<Integer> iterList = new ArrayList<>();
        ArrayList<Integer> positiveList = new ArrayList<>();
        ArrayList<Integer> falseList = new ArrayList<>();
        String nucleotideSequence = "";
        //Determine start and stop of sequence
        int Start = (int)(Math.random() * (d.size() - queryLength));
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
            } else {
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
            String[] out = BLAST.ProbabiliticBlast(s, i, d, wordSize, indexingThreshold, ungappedExtensionThreshold, delta, iterations, gapPenalty);
            {
                if (Double.valueOf(out[3]) >= Start && Double.valueOf(out[4]) <= Stop)
                {
                    positive += 1;
                }
                else
                {
                    falsePositive +=1;
                }
            }
            System.out.println("P: " + positive + " F: " + falsePositive);
            positiveList.add(positive);
            falseList.add(falsePositive);


            //falseList.add(falsePositive);
            iterList.add(testNumber);
            //Degrade sequence
            if (testNumber+1 < iterations)
            {
                String nucleotideSequence2 = "";
                for (int n =0;n<nucleotideSequence.length();n++)
                {
                    if (index.containsKey(n))
                    {
                        nucleotideSequence2 += (Math.random() < (d.getProbability(index.get(n))) - mutationAdjustement) ? String.valueOf(nucleotideSequence.charAt(n)) : nucleotide[rand.nextInt(4)];
                    }
                    else
                    {
                        nucleotideSequence2 += String.valueOf(nucleotideSequence.charAt(n));
                    }
                }
                nucleotideSequence = nucleotideSequence2;
                s = new Sequence(nucleotideSequence);
            }
        }
        Plot plt1 = Plot.create();
        Plot plt2 = Plot.create();
        plt1.plot().add(iterList, positiveList);
        plt2.plot().add(iterList,falseList);
        plt1.title("Positive matches");
        plt2.title("False Positive matches");
        plt1.xlabel("Mutation Rounds");
        plt2.xlabel("Mutation Rounds");
        plt1.ylabel("Found Sequences");
        plt2.ylabel("Found Sequences");
        //plt1.xticks(iterList);
        //plt2.xticks(iterList);
        try {
            plt1.show();
            plt2.show();
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
