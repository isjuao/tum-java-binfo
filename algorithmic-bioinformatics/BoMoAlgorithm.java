import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class BoMoAlgorithm {
    private ArrayList<String[]> readFiles(String name) {
        ArrayList<String[]> result = new ArrayList<>();
        try {
            File file = new File(name);
            BufferedReader reader = new BufferedReader(new FileReader(file));
            String s;
            while((s = reader.readLine())!= null) {
                String[] elements = s.split("\t");
                result.add(elements);
            }
            reader.close();
        } catch(IOException e) {
            System.out.println("An error occurred!");
            e.printStackTrace();
        }
        return result;
    }

    private char[] stringToChar(String s) {
        char[] c = new char[s.length()];
        for(int i = 0; i < s.length(); i++) {
            c[i] = s.charAt(i);
        }
        return c;
    }

    private ArrayList<Integer> algorithm(char[] t, int n, char[] s, int m) {
        int[] shift = new int[m + 1];
        shift = computeShiftTable(shift, s, m);
        ArrayList<Integer> iValues = new ArrayList<>();
        int i = 0;
        iValues.add(i);
        int j = m - 1;
        outer: while(i <= n - m) {
            while(t[i + j] == s[j]) {
                if(j == 0) {
                    // desired output format: m, n, i, i1-ik
                    iValues.add(0, iValues.get(iValues.size() - 1)); // insert 0 at start
                    iValues.add(0, n);  // insert n at start
                    iValues.add(0, m);  // insert m at start
                    return iValues;
                }
                j--;
            }
            i += shift[j];
            iValues.add(i);
            j = m - 1;
        }
        iValues.add(0, -1);
        iValues.add(0, n);
        iValues.add(0, m);
        return iValues;
    }

    private int[] computeShiftTable(int[] shift, char[] s, int m) {
        // Initialize
        for(int j = 0; j <= m; j++) {
            shift[j] = m;
        }
        // Part 1: sigma <= j
        int[] b = new int[m + 1];
        b[0] = -1;
        b[1] = 0;
        int i = 0;
        for(int l = 2; l <= m; l++) {
            while((i >= 0) && (s[m - i - 1] != s[m - l])) {
                int o = l - i - 1;
                shift[m - i - 1] = Math.min(shift[m - i - 1], o);
                i = b[i];
            }
            i++;
            b[l] = i;
        }
        // Part 2: sigma > j
        int j = 0;
        for(i = b[m]; i >= 0; i = b[i]) {
            int o = m - i;
            while(j < o) {
                shift[j] = Math.min(shift[j], o);
                j++;
            }
        }
        return shift;
    }

    public static void main(String[] args) {
        BoMoAlgorithm bm = new BoMoAlgorithm();
        // get all cases
        ArrayList<String[]> cases = bm.readFiles(args[0]);
        for(int i = 0; i < cases.size(); i++) {
            // convert to char arrays for easier use
            String[] strings = cases.get(i);
            char[] s = bm.stringToChar(strings[0]);
            char[] t = bm.stringToChar(strings[1]);
            ArrayList<Integer> values = bm.algorithm(t, t.length, s, s.length);
            String output = "";
            for(Integer value : values) {
                output += value + "\t";
            }
            System.out.println(output);
        }
    }
}
