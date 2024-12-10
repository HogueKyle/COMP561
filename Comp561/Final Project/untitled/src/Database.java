import java.util.Random;

public class Database
{
    private String database = "";

    /**
     * Create a random database
     * @param size Size of database to generate
     */
    public Database(int size)
    {
        Random rand = new Random();
        String[] nucleotide = {"A","G","C","T"};
        for (int i=0;i<size;i++)
        {
            database = database + nucleotide[rand.nextInt(4)];
        }
    }

    /**
     * Create a database from a string
     * @param sequence Sequence to make database around
     */
    public Database(String sequence)
    {
        database = sequence;
    }

    /**
     * Return a portion of the database
     * @param start Starting position of return sequence
     * @param end End position of return sequence
     * @return The subsequence string
     */
    public String subSequence(int start, int end)
    {
        return database.substring(start, end);
    }

    /**
     * Return the entire database
     * @return The database string
     */
    public String getDatabaseString()
    {
        return database;
    }

    /**
     * Return the size of the database
     * @return Size of the database
     */
    public int getDatabaseSize()
    {
        return database.length();
    }
    public void addRNucleotide(int size)
    {
        Random rand = new Random();
        String[] nucleotide = {"A","G","C","T"};
        for (int i=0;i<size;i++)
        {
            database = database + nucleotide[rand.nextInt(4)];
        }
    }
}
