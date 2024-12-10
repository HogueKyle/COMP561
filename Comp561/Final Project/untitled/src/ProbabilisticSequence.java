import java.util.ArrayList;

public class ProbabilisticSequence
{
    private ArrayList<ProbabilisticNucleotide> sequence = new ArrayList<ProbabilisticNucleotide>();
    /**
     * Create empty sequence
     */
    public ProbabilisticSequence(){}

    /**
     * Create sequence based off array list
     * @param sequence ArrayList of sequence
     */
    public ProbabilisticSequence(ArrayList<ProbabilisticNucleotide> sequence)
    {
        this.sequence = sequence;
    }
    private class ProbabilisticNucleotide
    {
        String nucleotide;
        double probability;
        int position;

        /**
         * Create a nucleotide object
         * @param nucleotide Base string
         * @param probability Probability double
         * @param position Position in sequence
         */

        public ProbabilisticNucleotide(String nucleotide, double probability, int position)
        {
            this.nucleotide = nucleotide;
            this.probability = probability;
            this.position = position;
        }
    }

    /**
     * Size of sequence
     * @return Size of sequence
     */
    public int size()
    {
        return sequence.size();
    }

    /**
     * Add to sequence
     * @param nucleotide Nucleotide to add
     * @param probability Probability to add
     * @param position Position of nucleotide to add
     */
    public void add(String nucleotide, double probability, int position)
    {
        sequence.add(new ProbabilisticNucleotide(nucleotide,probability,position));
    }

    /**
     * Return a subsequence of the nucleotide sequence
     * @param start Start position in sequence
     * @param stop End position in sequence
     * @return Subsequence of nucleotides
     */
    public ProbabilisticSequence subSequence(int start, int stop)
    {
        return new ProbabilisticSequence(new ArrayList<>(sequence.subList(start, stop)));
    }

    /**
     * Get a specific nucleotide in the sequence
     * @param position Position of nucleotide to return
     * @return String of nucleotide at specified position
     */
    public String getNucleotide(int position)
    {
        return sequence.get(position).nucleotide;
    }

    /**
     * Get the probability of a specific nucleotide in the sequence
     * @param position Position of nucleotide to return
     * @return Double of probability of nucleotide at specified position
     */
    public double getProbability(int position)
    {
        return sequence.get(position).probability;
    }
}
