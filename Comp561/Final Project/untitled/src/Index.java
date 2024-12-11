import java.util.ArrayList;
import java.util.Objects;

public class Index
{
    private ArrayList<ProbabilisticDatabase> words = new ArrayList<>();
    private int wordSize;

    /**
     * Index a sequence
     * @param d Sequence to index
     * @param wordSize Word size for each entry
     */
    public Index(ProbabilisticDatabase d, int wordSize)
    {
        this.wordSize = wordSize;
        for (int i = 0; i < d.size() - wordSize + 1; i++)
        {
            words.add(d.subDatabase(i, i+wordSize));
        }
    }

    /**
     * Check each word in index against a sequence
     * @param sequence Sequence to check
     * @param threshold Will show match if probabilistic score is above threshold
     * @return List of matching positions
     */
    public ArrayList<Integer> checkPosition(Sequence sequence, double threshold, int iterations)
    {
//        ArrayList<Integer> matches = new ArrayList<Integer>();
//        ArrayList<Sequence> randomSequence = new ArrayList<Sequence>(iterations);
//        for (int j=0;j<iterations;j++)
//        {
//            randomSequence.add(new Sequence(wordSize));
//        }
//        for(int x = 0; x < words.size(); x++)
//        {
//            if (BLAST.noisyScore(sequence, words.get(x), wordSize, iterations, randomSequence) >= threshold)
//            {
//                matches.add(x);
//            }
//        }
        ArrayList<Integer> matches = new ArrayList<Integer>();
        for(int x = 0; x < words.size(); x++)
        {
            if ((BLAST.calculateMatchScore(sequence, words.get(x), Boolean.TRUE, Double.MIN_VALUE)) >= (Math.log10(threshold)))
            {
                matches.add(x);
            }
        }
        return matches;
    }
    public String getWord(int x)
    {
        String r = "";
        for (int n=0;n<words.get(x).size();n++)
        {
            r += words.get(x).getNucleotide(n);
        }
        return r;
    }
}