import java.util.ArrayList;
import java.util.Arrays;

/**
 * Solves the three indistinguishable dice to two reddit post here:
 * https://www.reddit.com/r/math/comments/4eltj8/the_three_indistinguishable_dice_puzzle/
 *
 * Idea: 3 points on a circle with 6 points. Once you distinguish one of them, the distances clockwise to the
 * other two gives two independent die rolls. So, the method to distinguish is to first pick the one
 * with the minimum distance to its clockwise neighbor, and then use its value, and go clockwise that many times
 * The point you land on is the distinguished point.
 */
public class DiceRolls3To2
{
	public static int diff(int a, int b) {
		return (a - b + 6) % 6;
	}

	public static int[] order(int... arr) {
		int[] x = Arrays.copyOf(arr, arr.length);
		Arrays.sort(x);
		return x;
	}

	public static int get(int[] arr, int index) {
		return arr[index % arr.length];
	}

	public static boolean getBool(boolean[] arr, int index) {
		return arr[index % arr.length];
	}

	public static void main (String[] args) throws java.lang.Exception
	{
//		6 - 0, 0
//		30/2 - 0, x or x, x ->  5 * 2 + 5
//		20 -> a, b -> 5 * 4
//
//		aaa -> 0 0 = 6
		int[][] rolls = new int[6][6];
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				for (int k = 0; k < 6; k++) {
					int[] abc = order(i, j, k);
					int a = abc[0];
					int b = abc[1];
					int c = abc[2];
					int[] xyz = order(diff(b, a), diff(c, b), diff(a, c));

					int min = -100;

					boolean[] isMin = new boolean[3];
					for (int ind = 0; ind < 3; ind++) {
						isMin[ind] = diff(get(abc, ind + 1), get(abc, ind)) == xyz[0];
					}


					for (int ind = 0; ind < 3; ind++) {
						if (getBool(isMin,ind)) {
							if (!getBool(isMin, ind - 1 + 6) || getBool(isMin, ind + 1)) {
								min = ind;
							}
						}
					}
					min = (abc[min] + min) % 3;

					ArrayList<Integer> d = new ArrayList<>();
					d.add(diff(get(abc, min + 1), get(abc, min)));
					d.add(diff(get(abc, min + 2), get(abc, min)));

					int[] pair = order(d.get(0), d.get(1));
					rolls[pair[0]][pair[1]]++;
				}
			}
		}

		int total = 0;
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				if (i > j) {
					if (rolls[i][j] > 0) {
						System.out.println("error!!!!");
					} else {
						continue;
					}
				}
				total += rolls[i][j];
				System.out.println(i + "," + j + ":" + rolls[i][j]);
			}
		}
		System.out.println(total);
	}
}