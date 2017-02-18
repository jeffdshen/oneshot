import java.util.*;

/**
 * Problem: This is a generalization of DiceRolls3To2. Say probability distribution simulates another if we can map
 * events from one into the other such that we preserve the probailities. We then look at probability distributions
 * for simultaneous rolls of of k n-sided dice, which we denote kDn. In general, this question is NP-hard, but
 * we can provide two such mappings:
 *
 * Mapping 1: Takes k indistinguishable n-sided dice, and returns k-1 indistinguishable dice, so long as k divides n.
 *
 * Summary:
 * You plot each die modulo n (on a circle), and for each die, you take the clockwise distances with every other die
 * as a n-1-tuple. Choose the die with the smallest tuple, then if the value of the die is r, then choose the rth die
 * after that one, and take that die. Then your k-1 rolls are all other dice minus the latter die (mod n).
 *
 * Mapping 2: Takes k indistinguishable n-sided dice, and returns 2 indistinguishable dice, so long as k is even
 * and n is a power of 2.
 *
 * Summary: We first show that its true for n = 2. Then, we take a_i to be the number of dice with the value i.
 * Taking the sum of each half, we consider these as the number of heads and tails from kD2, and map these to 2D2.
 * The coin flip tells us the sum of first and second half of the mapped rolls (2Dn). In any of the 20, 02, or 11
 * cases, we can take the values and generate one die roll for each half, except for when one of the sides had a
 * sum of 0. In these cases, either they map to 20 or 02, in which case we can recurse, or they map to 11, in which
 * case we know theres a symmetry, so we can generate a distribution of 121 for each, and then combine them into 4
 * distinct groups.
 *
 * Lastly, we also try to make a program that finds a mapping given two probability distributions.
 *
 */
public class DiceRollsMapping {
    private static long pow(int n, int k) {
        long nk = 1;
        for (int i = 0; i < k; i++) {
            nk *= n;
        }
        return nk;
    }

    private static HashMap<ArrayList<Integer>, Integer> getDistribution(int n, int k) {
        long nk = pow(n, k);

        HashMap<ArrayList<Integer>, Integer> map = new HashMap<>();
        for (long i = 0; i < nk; i++) {
            ArrayList<Integer> dice = new ArrayList<>();
            long cur = i;
            for (int j = 0; j < k; j++) {
                dice.add((int) (cur % n));
                cur /= n;
            }

            Collections.sort(dice);
            if (!map.containsKey(dice)) {
                map.put(dice, 0);
            }

            // put in n, because the counts for k are n times the counts for k - 1.
            map.put(dice, map.get(dice) + 1);
        }

        return map;
    }

    private static void scale(Map<ArrayList<Integer>, Integer> dist, int factor) {
        for (Map.Entry<ArrayList<Integer>, Integer> entry : dist.entrySet()) {
            entry.setValue(entry.getValue() * factor);
        }
    }

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

    public static void mapDownBy1(int n, int k) {
        long nk = pow(n, k);

        HashMap<ArrayList<Integer>, Integer> count = new HashMap<>();

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
        Map<ArrayList<Integer>, Integer> expected = getDistribution(n, k - 1);
        scale(expected, n);
        System.out.println(expected.equals(count));
        System.out.println(count);
        System.out.println(expected);
    }

    static long calls = 0;

    private static boolean mapExistsHelper(
        ArrayList<Integer> a, ArrayList<Integer> b, int sum, int last, int index, HashSet<Integer> seen
    ) {
        calls++;
        if (calls % 10000 == 0) {
            System.out.println(calls + "," + sum + "," + index + "," + seen);
        }

        if (index >= b.size()) {
            return true;
        }

        int cur = -1;
        for (int i = 0; i < a.size(); i++) {
            if (seen.contains(i)) {
                continue;
            }

            if (cur >= 0 && Objects.equals(a.get(i), a.get(cur))) {
                continue;
            }

            cur = i;
            int nextLast = a.get(i);
            int next = sum + nextLast;
            if (next > b.get(index)) {
                continue;
            }

            if (last > 0 && nextLast > last) {
                continue;
            }

            seen.add(i);
            int nextIndex = index;
            if (next == b.get(index)) {
                next = 0;
                nextLast = 0;
                nextIndex++;
            }
            boolean exists = mapExistsHelper(a, b, next, nextLast, nextIndex, seen);
            seen.remove(i);

            if (exists) {
                return true;
            }
        }

        return false;
    }

    public static boolean mapExists(Collection<Integer> from, Collection<Integer> to) {
        ArrayList<Integer> a = new ArrayList<>(from);
        ArrayList<Integer> b = new ArrayList<>(to);
        Collections.sort(a);
        Collections.sort(b);
        Collections.reverse(a);
        System.out.println(a);
        System.out.println(b);
        calls = 0;
        return mapExistsHelper(a, b, 0, 0, 0, new HashSet<>());
    }

    public static boolean mapExists(int n, int k, int j) {
        HashMap<ArrayList<Integer>, Integer> from = getDistribution(n, k);
        HashMap<ArrayList<Integer>, Integer> to = getDistribution(n, j);
        scale(to, (int) pow(n, k - j));
        return mapExists(from.values(), to.values());
    }

    public static void main(String[] args) {
//        mapDownBy1(20, 5);
        System.out.println(mapExists(9, 4, 2));
    }

}
