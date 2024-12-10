import java.util.Random;

public class Sequence
{
    private StringBuilder sequence = new StringBuilder();

    /**
     * Create a random sequence
     * @param size Size of sequence to generate
     */
    public Sequence(int size)
    {
        Random rand = new Random();
        String[] nucleotide = {"A","G","C","T"};
        for (int i=0;i<size;i++)
        {
            sequence.append(nucleotide[rand.nextInt(4)]);
        }
    }

    /**
     * Create a sequence from a string
     * @param sequence Sequence to make sequence around
     */
    public Sequence(String sequence)
    {
        this.sequence.append(sequence);
    }
    public Sequence(int start, int stop, int size, String sequence)
    {
        Random rand = new Random();
        String[] nucleotide = {"A","G","C","T"};
        for (int i=0;i<start;i++)
        {
            sequence = sequence + nucleotide[rand.nextInt(4)];
        }
        sequence += sequence;
        for (int i=sequence.length();i<start;i++)
        {
            sequence = sequence + nucleotide[rand.nextInt(4)];
        }
    }

    /**
     * Return a portion of the sequence
     * @param start Starting position of return sequence
     * @param end End position of return sequence
     * @return The subsequence string
     */
    public Sequence subSequence(int start, int end)
    {
        return new Sequence(sequence.substring(start, end));
    }

    /**
     * Return the entire sequence
     * @return The sequence string
     */
    public String getNucleotide(int p)
    {
        return String.valueOf(sequence.charAt(p));
    }
    public Sequence addZero()
    {
        return new Sequence(sequence.insert(0, '0').toString());
    }
//    public String getSequenceString()
//    {
//        return sequence;
//    }

    /**
     * Return the size of the sequence
     * @return Size of the sequence
     */
    public int size()
    {
        return sequence.length();
    }
    public void addRNucleotide(int size)
    {
        Random rand = new Random();
        String[] nucleotide = {"A","G","C","T"};
        for (int i=0;i<size;i++)
        {
            sequence.append(nucleotide[rand.nextInt(4)]);
        }
    }
}
