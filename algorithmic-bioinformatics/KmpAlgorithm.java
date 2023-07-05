import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class KmpAlgorithm {
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

    private float[] algorithm(char[] t, int n, char[] s, int m) {
        long startTime = System.nanoTime();
        // initialize
        int[] b = new int[m + 1];
        this.borders(b, m, s);
        int i = 0;
        int j = 0;
        int x = -1;
        // step through text, check if search string would still fit
        outer: while(i <= n - m) {
            // check if match
            while(t[i + j] == s[j]) {
                j++;
                // check if checked whole search string
                if(j == m) {
                    // found search string
                    x = i;
                    break outer;
                }
            }
            // mismatch: shift search string to the right in text
            i = i + (j - b[j]);
            // mismatch: start comparing after boarder (already matched)
            j = Math.max(0, b[j]);
        }
        // finish
        int test;
        long endTime = System.nanoTime();
        float[] result = new float[4];
        result[0] = m;
        result[1] = n;
        result[2] = x;
        result[3] = (float)((endTime - startTime) / 1_000_000_000.0);
        return result;
    }

    private void borders(int[] b, int m, char[] s) {
        b[0] = -1;
        b[1] = 0;
        int i = 0;
        for(int j = 2; j <= m; j++) {
            // check if there actually is (not saying we can, but just the possibility)
            // (1) "eigentlicher Rand vom letzten Rand (blau)" to expand (I think) and
            // (2) we cannot expand "bekannter eigentlicher Rand (grün)"
            while((i >= 0) && (s[i] != s[j-1])) {
                // place i after/at "eigentlicher Rand vom letzten Rand" (blau) to check if we really can expand
                i = b[i];
            }
            // we match and can expand "bekannter eigentlicher Rand (grün)"
            i++;
            // found our border and can move on in search string
            b[j] = i;
        }
    }

    public static void main(String[] args) {
        KmpAlgorithm kmp = new KmpAlgorithm();
        // get all cases
        ArrayList<String[]> cases = kmp.readFiles(args[0]);
        for(int i = 0; i < cases.size(); i++) {
            // convert to char arrays for easier use
            String[] strings = cases.get(i);
            char[] s = kmp.stringToChar(strings[0]);
            char[] t = kmp.stringToChar(strings[1]);
            float[] values = kmp.algorithm(t, t.length, s, s.length);
            System.out.println(
                    (int) values[0] + "\t" + (int) values[1] + "\t" + (int) values[2] + "\t" + values[3]
            );
        }
    }
}
