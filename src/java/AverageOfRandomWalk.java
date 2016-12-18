import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.*;

/**
 * Problem: You take a random walk on a number line starting at the origin for n steps,
 * moving left or right by one with equal probability.
 *
 * Given the average position, what is the expected final position?
 *
 * Note: the reverse is easy - if the final position is x, then the expected average position is just x/2.
 *
 * Method: Let t be the number of steps. Use sum of positions instead of average (easier to work with).
 * First step is either up or down, and that affects the sum of positions by either t or -t.
 * We can now recurse to t-1 steps and the new sum (either adjusted by -t or t respectively).
 * Finding the sum and count of final positions for the two t-1 cases, we can figure out the sum and count for t.
 * The recursion can be flipped around into an iterative dp solution.
 *
 * Runtime: O(t^3), t can go up to around 400 without taking too long.
 *
 * Conclusion: The ratio of the expected final position and x starts as 1.5 for small x and then goes to 2.
 * When plotted against each other, the curve is convex, and looks almost like a half circle (but not exactly).
 *
 * Uncomment out the commented code in main to get the points in tabular format.
 */
public class AverageOfRandomWalk {
    public static class Avg {
        public final BigInteger sum;
        public final BigInteger count;

        public Avg(BigInteger sum, BigInteger count) {
            this.sum = sum;
            this.count = count;
        }

        public double getAvg() {
            return new BigDecimal(sum).divide(new BigDecimal(count), 20, RoundingMode.HALF_EVEN).doubleValue();
        }

        @Override
        public String toString() {
            return "Avg{" +
                "sum=" + sum +
                ", count=" + count +
                '}';
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Avg avg = (Avg) o;

            if (sum != null ? !sum.equals(avg.sum) : avg.sum != null) return false;
            return count != null ? count.equals(avg.count) : avg.count == null;

        }

        @Override
        public int hashCode() {
            int result = sum != null ? sum.hashCode() : 0;
            result = 31 * result + (count != null ? count.hashCode() : 0);
            return result;
        }
    }

    private static BigInteger getSum(BigInteger[] arr, int sum, int t) {
        if (sum < 0) {
            return getSum(arr, -sum, t).negate();
        }

        if (sum > t * (t + 1) / 2) {
            return BigInteger.valueOf(0);
        }

        return arr[sum];
    }

    private static BigInteger getCount(BigInteger[] arr, int sum, int t) {
        if (sum < 0) {
            return getCount(arr, -sum, t);
        }

        if (sum > t * (t + 1) / 2) {
            return BigInteger.ZERO;
        }

        return arr[sum];
    }

    /**
     * Use iterative dp to return the sum of positions (divide by t or t + 1 to get avg) -> expected final position
     * @param t how many steps to take
     * @return a map from sum of positions to expected final position
     */
    public static Map<Integer, Avg> dp(int t) {
        int max = t * (t + 1) / 2;
        BigInteger[] S = new BigInteger[max + 1]; // sums
        BigInteger[] C = new BigInteger[max + 1]; // counts

        BigInteger[] nextS = new BigInteger[max + 1];
        BigInteger[] nextC = new BigInteger[max + 1];

        for (int i = 0; i <= max; i++) {
            S[i] = BigInteger.ZERO;
            C[i] = BigInteger.ZERO;
            nextS[i] = BigInteger.ZERO;
            nextC[i] = BigInteger.ZERO;
        }

        S[1] = BigInteger.ONE;
        C[1] = BigInteger.ONE;
        for (int i = 2; i <= t; i++) {
            for (int j = 0; j <= i * (i + 1) / 2; j++) {
                int down = j + i;
                int up = j - i;
                BigInteger downS = getSum(S, down, t - 1);
                BigInteger downC = getCount(C, down, t - 1);

                BigInteger upS = getSum(S, up, t - 1);
                BigInteger upC = getCount(C, up, t - 1);

                nextS[j] = downS.subtract(downC).add(upS).add(upC);
                nextC[j] = upC.add(downC);
            }
            BigInteger[] tempS = S;
            S = nextS;
            nextS = tempS;

            BigInteger[] tempC = C;
            C = nextC;
            nextC = tempC;
        }

        HashMap<Integer, Avg> ans = new HashMap<>();
        for (int i = 0; i <= max; i++) {
            ans.put(i, new Avg(S[i], C[i]));
            ans.put(-i, new Avg(S[i].negate(), C[i]));
        }

        return ans;
    }

    public static void main(String[] args) {
        int n = 110;

        Map<Integer, Avg> ans = dp(n);
        double min = 100; // should be at most 2
        for (int i = -n * (n + 1) / 2; i <= n * (n + 1) / 2; i++) {
            Avg fin = ans.get(i);

            if (fin.count.equals(BigInteger.ZERO)) {
                continue;
            }

            double avg = i * 1.0 / (n + 1);
            double ratio = avg == 0d ? min : fin.getAvg() / avg;
            System.out.println(
                String.format("%f\t%f",
                    avg,
                    fin.getAvg()
                )
            );
            min = Math.min(min, ratio);
        }
//        System.out.println(min);
    }
}
