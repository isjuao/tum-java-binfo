import java.util.HashMap;
import java.util.Hashtable;

public class Ebc {
    public HashMap alphabet;

    private int[][] calcEBC(char[] s, int m) {
        this.alphabet = new HashMap<>();
        for(int i = 0; i < m; i++) {
            if(!alphabet.containsKey(s[i])) {
                alphabet.put(s[i], 0);
            }
        }
        int[][] ebcTable = new int[alphabet.size()][m];
        int counter = 0;
        for(Object a : alphabet.keySet()) {
            ebcTable[counter][0] = -1;
            for(int k = 1; k < m; k++) {
                if(s[k - 1] == (char)a) {
                    ebcTable[counter][k] = k - 1;
                } else {
                    ebcTable[counter][k] = ebcTable[counter][k - 1];
                }
            }
            counter++;
        }
        return ebcTable;
    }

    public static void main(String[] args) {
        Ebc ebc = new Ebc();
        char[] s = {'a', 'b', 'a', 'c', 'c', 'b', 'a', 'a'};
        int[][] ebcTable = ebc.calcEBC(s, s.length);
        int i = 0;
        for(Object a : ebc.alphabet.keySet()) {
            System.out.print(a + ": ");
            for(int j = 0; j < ebcTable[i].length; j++) {
                System.out.print(ebcTable[i][j] + "\t");
            }
            System.out.print("\n");
            i++;
        }
    }
}
