import java.awt.*;
import java.awt.event.KeyEvent;

/**
 * Simple script for keyboard macro cause apparently nothing else can do key down, hold for a bit, then key up
 */
public class HyperLightDrifterDash {
    static long last = 0;

    public static void press(Robot robot, double hold, double delay) {
        robot.keyPress(KeyEvent.VK_SPACE);
        robot.delay((int)(hold / 1000 + 1));
        robot.keyRelease(KeyEvent.VK_SPACE);
        robot.delay((int)(delay / 1000 + 1));
        long temp = System.currentTimeMillis();
        System.out.println(temp - last);
        last = temp;
    }

    public static void main(String[] args) throws AWTException {
        Robot robot = new Robot();
        robot.delay(3000);
        int frame = 16666;
        double after = 160000 / frame;
        double hold = 100000;
        press(robot, hold, frame * 14);
        press(robot, hold, frame * 13);
        press(robot, hold, frame * 12);
        press(robot, hold, frame * 11);
        press(robot, hold, frame * 10);

        press(robot, hold, frame * after);
        press(robot, hold, frame * after);
        press(robot, hold, frame * after);

        for (int i = 0; i < 800; i++) {
            press(robot, hold, frame * after);
        }
    }
}
