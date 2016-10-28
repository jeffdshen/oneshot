import java.util.*;

/**
 * https://www.research.ibm.com/haifa/ponderthis/challenges/October2016.html
 */
public class PonderThisSep2016 {
    public static ArrayList<Integer> stack = new ArrayList<>();
    public static ArrayList<Integer> pool = getPool();
    public static Set<Integer> seen = new HashSet<>();
    public static Map<Integer, Integer> viol = getViolations();

    public static ArrayList<Integer> getPool() {
        ArrayList<Integer> pool = new ArrayList<>();
        for (int i = 0; i < 10000; i++) {
            if (check(i)) {
                pool.add(i);
            }
        }

        return pool;
    }

    public static boolean check(int num) {
        boolean[] seen = new boolean[5];
        for (int i = 0; i < 4; i++) {
            if (num % 10 >= 5 || num % 10 == 0) {
                return false;
            }

            if (seen[num % 10]) {
                return false;
            }

            seen[num % 10] = true;
            num /= 10;
        }

        return true;
    }

    public static int violations() {
        int x = stack.size() - 1;
        int y = stack.size() - 2;
        int z = stack.size() - 3;
        int count = 0;
        if (y >= 0) {
            count += violation(stack.get(x), stack.get(y));
        }

        if (z >= 0) {
            count += violation(stack.get(x), stack.get(z));
        }

        return count;
    }

    public static int violationHelper(int x, int y) {
        int count = 0;
        for (int i = 0; i < 4; i++) {
            if (x % 10 == y % 10) {
                count++;
            }
            x /= 10;
            y /= 10;
        }

        return count;
    }

    public static Map<Integer,Integer> getViolations() {
        Map<Integer, Integer> viol = new HashMap<>();
        for (int i : pool) {
            for (int j : pool) {
                viol.put(i * 10000 + j, violationHelper(i, j));
            }
        }

        return viol;
    }

    public static int violation(int x, int y) {
        return viol.get(x * 10000 + y);
    }

    public static boolean dfs(int depth, int v, int max) {
        if (v > max) {
            return false;
        }

        if (depth == 24) {
            System.out.println(v);
            for (int x : stack) {
                System.out.println(x);
            }
            return true;
        }

        for (int x : pool) {
            if (seen.contains(x)) {
                continue;
            }

            seen.add(x);
            stack.add(x);
            boolean works = dfs(depth + 1, v + violations(), max);
            seen.remove(x);
            stack.remove(depth);
            if (works) {
                return true;
            }
        }

        return false;
    }

    public static void main(String[] args) {
        for (int i = 0; i < 10000; i++) {
            boolean works = dfs(0, 0, i);
            if (works) {
                break;
            }
        }
    }
}
