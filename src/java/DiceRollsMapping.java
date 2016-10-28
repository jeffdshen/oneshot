import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

/**
 * Takes k indistinguishable n-sided dice, and returns k-1 indistinguishable dice. Precondition: k divides n.
 *
 * Summary:
 * You plot each die modulo n (on a circle), and for each die, you take the clockwise distances with every other die
 * as a n-1-tuple. Choose the die with the smallest tuple, then if the value of the die is r, then choose the rth die
 * after that one, and take that die. Then your k-1 rolls are all other dice minus the latter die (mod n).
 *
 */
public class DiceRollsMapping {
    private static int diff(int a, int b, int n) {
        return (a - b + n) % n;
    }

    private static int compareTo(ArrayList<Integer> a, ArrayList<Integer> b) {
        if (a.size() != b.size()) {
            throw new IllegalArgumentException();
        }

        for (int i = 0; i < a.size(); i++) {
            if (!a.get(i).equals(b.get(i))) {
                return a.get(i) - b.get(i);
            }
        }

        return 0;
    }

    /**
     * Takes a sorted arraylist with k rolls, and returns a sorted arraylist with k - 1 rolls.
     * All numbers are from 0 to n - 1.
     */
    private static ArrayList<Integer> roll(ArrayList<Integer> dice, int n) {
        ArrayList<ArrayList<Integer>> tuples = new ArrayList<>();

        int min = 0;
        for (int i = 0; i < dice.size(); i++) {
            ArrayList<Integer> tuple = new ArrayList<>();
            for (int j = 1; j < dice.size(); j++) {
                tuple.add(diff(dice.get((j + i) % dice.size()), dice.get(i), n));
            }

            Collections.sort(tuple); // just in case, but this should already be sorted.
            tuples.add(tuple);
            if (compareTo(tuple, tuples.get(min)) < 0) {
                min = i;
            }
        }

        return tuples.get((min + dice.get(min)) % dice.size());
    }

    public static void main(String[] args) {
        int n = 20;
        int k = 5;

        long nk = 1;
        for (int i = 0; i < k; i++) {
            nk *= n;
        }

        long nk1 = 1;
        for (int i = 0; i < k - 1; i++) {
            nk1 *= n;
        }

        HashMap<ArrayList<Integer>, Integer> count = new HashMap<>();
        HashMap<ArrayList<Integer>, Integer> expected = new HashMap<>();

        for (long i = 0; i < nk; i++) {
            ArrayList<Integer> dice = new ArrayList<>();
            long cur = i;
            for (int j = 0; j < k; j++) {
                dice.add((int) (cur % n));
                cur /= n;
            }

            Collections.sort(dice);
            ArrayList<Integer> roll = roll(dice, n);
            if (!count.containsKey(roll)) {
                count.put(roll, 0);
            }

            count.put(roll, count.get(roll) + 1);
        }

        // verify
        for (long i = 0; i < nk1; i++) {
            ArrayList<Integer> dice = new ArrayList<>();
            long cur = i;
            for (int j = 0; j < k - 1; j++) {
                dice.add((int) (cur % n));
                cur /= n;
            }

            Collections.sort(dice);
            if (!expected.containsKey(dice)) {
                expected.put(dice, 0);
            }

            // put in n, because the counts for k are n times the counts for k - 1.
            expected.put(dice, expected.get(dice) + n);
        }
        System.out.println(expected.equals(count));
        System.out.println(count);
        System.out.println(expected);
    }

}
