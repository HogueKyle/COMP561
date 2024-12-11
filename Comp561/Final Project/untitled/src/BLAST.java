import java.util.ArrayList;
import java.util.stream.IntStream;

public class BLAST
{
    public static double safeMultiply(double a, double b, double c)
    {
        if ((a == 0) && (b == 0))
        {
            return Math.log10(c);
        } else if (b==0)
        {
            return a + Math.log10(c);
        }
        else
        {
            return a + Math.log10(b);
        }
    }
    public static double calculateMatchScore(Sequence s, ProbabilisticDatabase d, Boolean operation, double safetyZero)
    {
        double score;
        if (operation)
        {
            score = 0; // 1;
        }
        else
        {
            score= 0;
        }
        for(int c=0;c<s.size();c++)
        {
            if (s.getNucleotide(c).equals(d.getNucleotide(c)))
            {
                if (operation)
                {
                    score = safeMultiply(score, d.getProbability(c), safetyZero);
                }
                else
                {
                    score += d.getProbability(c);
                }
            }
            else
            {
                if (operation)
                {
                    score = safeMultiply(score, ((1 - d.getProbability(c)) / 3), safetyZero);
                }
                else
                {
                    score += (1 - d.getProbability(c)) / 3;
                }
            }
        }
        return score;
    }

    public static double noisyScore(Sequence s, ProbabilisticDatabase d, int size, int iterations, ArrayList<Sequence> randomSequences, Boolean operation, double safetyZero)
    {

        double noise = 0;
        //Correct random sequence size
        if (randomSequences.get(0).size() > size)
        {
            for (int i = 0; i < iterations; i++)
            {
                noise += calculateMatchScore(randomSequences.get(i).subSequence(0,size), d, operation, safetyZero);
            }
        } else if (randomSequences.get(0).size() < size)
        {
            for (int i = 0; i < iterations; i++)
            {
                randomSequences.get(i).addRNucleotide(size - randomSequences.get(i).size());
                noise += calculateMatchScore(randomSequences.get(i), d, operation, safetyZero);
            }
        }
        else
        {
            for (int i = 0; i < iterations; i++)
            {
                noise += calculateMatchScore(randomSequences.get(i), d, operation, safetyZero);
            }
        }
        noise /= iterations;
        return calculateMatchScore(s, d, operation, safetyZero)-noise;
        //return Math.pow(10,calculateMatchScore(s, d, operation, safetyZero))/Math.pow(10,noise);
    }

    /**
     * Perform BLAST for a probabilistic sequence
     * @param s Sequence
     * @param i Index
     * @param d Database
     * @param wordSize Word size in index
     * @param indexingThreshold Threshold to match index to sequence
     * @param ungappedExtensionThreshold Threshold to proceed from ungapped extension to gapped extension
     * @param delta Stopping condition for ungapped extension
     * @param iterations Number of noise sequences to create to find a mean
     */
    public static String[] ProbabiliticBlast(Sequence s, Index i, ProbabilisticDatabase d, int wordSize, double indexingThreshold, double ungappedExtensionThreshold, double delta, int iterations, double gapPenalty, double safetyZero)
    {
        int matchSize = 0;
        ArrayList<String[]> NWMatch = new ArrayList<>();
        //Phase 1
        for(int sequencePosition = 0; sequencePosition < s.size() - wordSize + 1; sequencePosition++)
        {
            //System.out.println("Position " + sequencePosition + " ----------------------------------------------------------------------------------------------------------------");
            //Check for matches against database index
            ArrayList<Integer> matches = i.checkPosition(s.subSequence(sequencePosition, sequencePosition + wordSize), indexingThreshold, iterations);
            matchSize += matches.size();
            //If matches proceed to ungapped extension
            for (int m = 0; m < matches.size(); m++)
            {
                ArrayList<Sequence> randomSequence = new ArrayList<Sequence>(iterations);
                for (int j=0;j<iterations;j++)
                {
                    randomSequence.add(new Sequence(i.getWord(matches.get(m))));
                }
                int StartDF = matches.get(m);
                int StopDF = StartDF + wordSize;
                int StartSF = sequencePosition;
                int StopSF = sequencePosition + wordSize;
                //System.out.println("Index Match Found \uD83D\uDE0E! Start Sequence: " + StartSF + " Stop Sequence: " + StopSF + " Start Database: " + StartDF + " Stop Database: " + StopDF);
                //System.out.println("Proceeding to ungapped extension \uD83D\uDE00:");
                int bestStartS = StartSF;
                int bestStopS = StopSF;
                int bestStartD = StartDF;
                int bestStopD = StopDF;
                double prob = 0;
                //Initial score
                double initialScore = noisyScore(s.subSequence(StartSF, StopSF), d.subDatabase(StartDF, StopDF), StopDF - StartDF, iterations, randomSequence, Boolean.TRUE, safetyZero);
                double bestScore = initialScore;
                //Forward extension
                do
                {
                    StopSF += 1;
                    StopDF += 1;
                    if (StopSF > s.size() || StopDF > d.size())
                    {
                        StopSF--;
                        StopDF--;
                        break;
                    }
                    prob = noisyScore(s.subSequence(StartSF, StopSF), d.subDatabase(StartDF, StopDF),StopDF - StartDF, iterations, randomSequence, Boolean.TRUE, safetyZero);
                    if (prob > bestScore)
                    {
                        bestScore = prob;
                        bestStopD = StopDF;
                        bestStopS = StopSF;
                    }
                }while (prob >= bestScore - delta);
                //Backwards extension
                int StartDB = matches.get(m);
                int StopDB = StartDF + wordSize;
                int StartSB = sequencePosition;
                int StopSB = sequencePosition + wordSize;
                bestScore = initialScore;
                do
                {
                    StartSB -= 1;
                    StartDB -= 1;
                    if (StartSB < 0 || StartDB < 0)
                    {
                        StartSB++;
                        StartDB++;
                        break;
                    }
                    prob = noisyScore(s.subSequence(StartSB, StopSB), d.subDatabase(StartDB, StopDB),StopDB - StartDB, iterations, randomSequence, Boolean.TRUE, safetyZero);
                    if (prob > bestScore)
                    {
                        bestScore = prob;
                        bestStartS = StartSB;
                        bestStartD = StartDB;
                    }
                }while (prob >= bestScore - delta/*prob >= delta*/);
                //Final best score
                prob = noisyScore(s.subSequence(bestStartS, bestStopS), d.subDatabase(bestStartD, bestStopD),bestStopD - bestStartD, iterations, randomSequence, Boolean.TRUE, safetyZero);
                //prob = calculateMatchScore(s.subSequence(bestStartS, bestStopS), d.subDatabase(bestStartD, bestStopD), Boolean.TRUE, safetyZero);
                if (prob >= ungappedExtensionThreshold)
                {
                    System.out.println("Ungapped Extension Match Found \uD83E\uDDEC! Start Sequence: " + bestStartS + " Stop Sequence: " + bestStopS + " Start Database: " + bestStartD + " Stop Database: " + bestStopD + " Score: " + (prob));
                    //NW
                    System.out.println("NW " + StartDB + " to " + StopDF);
                    NWMatch.add(NeedlemanWunsch.NW(s.subSequence(0,s.size()), d.subDatabase(matches.get(m) - sequencePosition, matches.get(m) - sequencePosition + s.size()), StartDB, StopDF, gapPenalty));
                }
                else
                {
                    //System.out.println("Ungapped Extension Failed \uD83D\uDE11: " + prob);
                }
            }
        }
        double bestScore = Double.MIN_VALUE;
        int bestIndex = Integer.MIN_VALUE;
        //Return best scoring NW
        for(int z=0;z<NWMatch.size();z++)
        {
            if (Double.valueOf(NWMatch.get(z)[2]) > bestScore)
            {
                bestScore = Double.valueOf(NWMatch.get(z)[2]);
                bestIndex = z;
            }
        }
        //Report results
        if (bestIndex != Integer.MIN_VALUE) {
            String[] re = NWMatch.get(bestIndex);
            re[5] = Integer.toString(matchSize);
            re[6] = Integer.toString(NWMatch.size());
            return re;
        }
        else
        {
            return new String[]{"0", "0", "0", "0", "0", "0", "0"};
        }
    }
}
