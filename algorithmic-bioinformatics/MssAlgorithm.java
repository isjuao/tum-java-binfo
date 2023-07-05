import java.io.*;
import java.util.ArrayList;

public class MSS_clever {
    private ArrayList<Integer[]> readFiles(String name) {
        ArrayList<Integer[]> result = new ArrayList<>();
        try {
            File file = new File(name);
            BufferedReader reader = new BufferedReader(new FileReader(file));
            String s;
            while((s = reader.readLine())!= null) {
                String[] elements = s.split("\t");
                Integer[] a = new Integer[elements.length];
                for(int i = 0; i < elements.length ; i++) {
                    a[i] = Integer.parseInt(elements[i]);
                }
                result.add(a);
            }
            reader.close();
        } catch(IOException e) {
            System.out.println("An error occurred!");
            e.printStackTrace();
        }
        return result;
    }

    private float[] algorithm(int[] a, int n) {
        long startTime = System.nanoTime();
        int maxscore = 0;
        int rmaxscore = 0;
        int l = 1;
        int r = 0;
        int rstart = 1;
        for(int i = 1; i <= n; i++) {
            if(rmaxscore > 0) {
                rmaxscore = rmaxscore + a[i-1];
            } else {
                rmaxscore = a[i-1];
                rstart = i;
            }
            if(rmaxscore > maxscore) {
                maxscore = rmaxscore;
                l = rstart;
                r = i;
            }
        }
        long endTime = System.nanoTime();
        float[] result = new float[5];
        result[0] = a.length;
        result[1] = maxscore;
        result[2] = l;
        result[3] = r;
        result[4] = (float)((endTime - startTime) / 1_000_000_000.0);
        return result;
    }

    public static void main(String[] args) {
        MSS_clever mss = new MSS_clever();
        // int[] a = {5,-2,5,-2,1,-9,5,-2,4,-5,1,-2,3,-1,5,-3,2,-1,2};
        String name;
        ArrayList<Integer[]> cases = mss.readFiles(args[0]);
        for(int i = 0; i < cases.size(); i++) {
            Integer[] aInt = cases.get(i);
            int[] a = new int[aInt.length];
            for(int j = 0; j < aInt.length; j++) {
                a[j] =  aInt[j];
            }
            float[] values = mss.algorithm(a, a.length);
            System.out.println((int)values[0] + "\t" + (int)values[1] + "\t" + (int)values[2] + "\t" + (int)values[3]
                    + "\t" + values[4]);
        }
    }
}
