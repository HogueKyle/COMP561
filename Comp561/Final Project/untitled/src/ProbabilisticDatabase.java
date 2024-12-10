import java.util.ArrayList;

public class ProbabilisticDatabase
{
    private ArrayList<ProbabilisticNucleotide> database = new ArrayList<ProbabilisticNucleotide>();
    /**
     * Create empty database
     */
    public ProbabilisticDatabase(){}

    /**
     * Create database based off array list
     * @param database ArrayList of database
     */
    public ProbabilisticDatabase(ArrayList<ProbabilisticNucleotide> database)
    {
        this.database = database;
    }
    private class ProbabilisticNucleotide
    {
        String nucleotide;
        double probability;

        /**
         * Create a nucleotide object
         * @param nucleotide Base string
         * @param probability Probability double
         */

        public ProbabilisticNucleotide(String nucleotide, double probability)
        {
            this.nucleotide = nucleotide;
            this.probability = probability;
        }
    }


    /**
     * Size of database
     * @return Size of database
     */
    public int size()
    {
        return database.size();
    }

    /**
     * Add to database
     * @param nucleotide Nucleotide to add
     * @param probability Probability to add
     */
    public void add(String nucleotide, double probability)
    {
        database.add(new ProbabilisticNucleotide(nucleotide,probability));
    }

    /**
     * Return a subdatabase of the nucleotide database
     * @param start Start position in database
     * @param stop End position in database
     * @return Subdatabase of nucleotides
     */
    public ProbabilisticDatabase subDatabase(int start, int stop)
    {
        return new ProbabilisticDatabase(new ArrayList<>(database.subList(start, stop)));
    }
    public ProbabilisticDatabase addZero()
    {
        ArrayList<ProbabilisticNucleotide> database2 = (ArrayList<ProbabilisticNucleotide>) database.clone();
        database2.add(0, new ProbabilisticNucleotide("0",0));
        return new ProbabilisticDatabase(database2);
    }
    /**
     * Get a specific nucleotide in the database
     * @param position Position of nucleotide to return
     * @return String of nucleotide at specified position
     */
    public String getNucleotide(int position)
    {
        return database.get(position).nucleotide;
    }

    /**
     * Get the probability of a specific nucleotide in the database
     * @param position Position of nucleotide to return
     * @return Double of probability of nucleotide at specified position
     */
    public double getProbability(int position)
    {
        return database.get(position).probability;
    }

}
