//TIP To <b>Run</b> code, press <shortcut actionId="Run"/> or
// click the <icon src="AllIcons.Actions.Execute"/> icon in the gutter.

public class NeedlemanWunsch
{
    public static String[] NW(Sequence s, ProbabilisticDatabase t, int startIndex, int stopIndex, double gap)
    {
        double up;
        double left;
        double diagonal;

        double[][] X = new double[s.size() + 1][t.size() + 1];
        String[][] arrow = new String[s.size() + 1][t.size() + 1];
        //APPEND 0 to STRING
        t = t.addZero();
        s = s.addZero();
        for (int i = 0; i < (s.size()); i++)
        {
            for (int j = 0; j < (t.size()); j++)
            {
                up = Integer.MIN_VALUE;
                left = Integer.MIN_VALUE;
                diagonal = Integer.MIN_VALUE;
                //COMPARE CASES
                //Zero cases
                if (t.getNucleotide(j).equals("0") && s.getNucleotide(i).equals("0"))
                {
                    X[i][j] = 0;
                } else
                {
                    //Up
                    if (i - 1 >= 0)
                    {
                        up = X[i - 1][j] + gap;
                    }
                    //Left
                    if (j - 1 >= 0)
                    {
                        left = X[i][j - 1] + gap;
                    }
                    //Diagonal
                    if ((i - 1 >= 0) && (j - 1 >= 0))
                    {
                        diagonal = X[i - 1][j - 1] + ((s.getNucleotide(i).equals(t.getNucleotide(j))) ? t.getProbability(j) : (1 - t.getProbability(j))/3);
                                //M[MIndex(s.charAt(i))][MIndex(t.charAt(j))];
                    }
                    X[i][j] = Double.max(Double.max(up, left), diagonal);
                    if ((X[i][j]) == up)
                    {
                        arrow[i][j] = "up";
                    }
                    if ((X[i][j]) == left)
                    {
                        arrow[i][j] = "left";
                    }
                    if ((X[i][j]) == diagonal)
                    {
                        arrow[i][j] = "diagonal";
                    }
                }
            }
        }
        //Retrieve optimal allignment
        int i = s.size() -1;
        int j = t.size() -1;
        String alignmentS = " ";
        String alignmentT = " ";
        int firstDiagonal = 0;
        int lastDiagonal = 0;
        Boolean diagonalFlag = Boolean.TRUE;
        while (i != 0 || j != 0)
        {
            switch (arrow[i][j])
            {
                case "diagonal":
                    alignmentS = s.getNucleotide(i) + alignmentS;
                    alignmentT = t.getNucleotide(j) + alignmentT;
                    i -= 1;
                    j -= 1;
                    if (diagonalFlag)
                    {
                        diagonalFlag = Boolean.FALSE;
                        lastDiagonal = j + startIndex + 1;
                    }
                    firstDiagonal = j + startIndex;
                    break;
                case "left":
                    alignmentS = "-" + alignmentS;
                    alignmentT = t.getNucleotide(j) + alignmentT;
                    j -= 1;
                    break;
                case "up":
                    alignmentS = s.getNucleotide(i) + alignmentS;
                    alignmentT = "-" + alignmentT;
                    i -= 1;
                    break;
            }
        }
        alignmentS = alignmentS.replaceAll("\\s+","");
        alignmentT = alignmentT.replaceAll("\\s+","");
        System.out.println(alignmentS);
        System.out.println(alignmentT);
        System.out.println(X[s.size()-1][t.size()-1]);
        String[] r = {Integer.toString(firstDiagonal), Integer.toString(lastDiagonal), Double.toString(X[s.size()-1][t.size()-1]), alignmentS, alignmentT, "", ""};
        return r;
    }
}