import java.io.File;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.lang.reflect.Array;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;
public class Main {



    public static void main(String[] args) {
        //Read sequence and probability information
        //Read file
        try {
            URL databaseURL = new URL("https://www.cs.mcgill.ca/~blanchem/561/probabilisticGenome/chr22.maf.ancestors.42000000.complete.boreo.fa");
            URL probabilityURL = new URL("https://www.cs.mcgill.ca/~blanchem/561/probabilisticGenome/chr22.maf.ancestors.42000000.complete.boreo.conf");
            Scanner databaseReader = new Scanner(databaseURL.openStream());
            Scanner probabilityReader = new Scanner(probabilityURL.openStream());
            String probabilitiesUnbound = databaseReader.next();
            ProbabilisticDatabase d = new ProbabilisticDatabase();
            //Create sequence of nucleotide, probability and object
            int position = 0;
            while (probabilityReader.hasNext()) {
                d.add(String.valueOf(probabilitiesUnbound.charAt(position)), probabilityReader.nextDouble());
                position++;
            }
            databaseReader.close();
            probabilityReader.close();
            //Test.Run(3,0.05,100,d,0,10,10,3.5,80,1,-1);
            Test.Run(1,0.0,100,d,0,10,10,0.1,80,1,-1);
            //Same sequence

            //BLAST.ProbabiliticBlast(seq, i, d, 11, 0.5, 0.5, 5, 10);

        } catch (FileNotFoundException e) {
            System.out.println("Could not find file \uD83D\uDE31");
        } catch (IOException e) {
            System.out.println("Could not find file at URL \uD83D\uDE31");
        }
    }
}