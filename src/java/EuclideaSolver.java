import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.*;

/**
 * DFS for compass and straight-edge drawings.
 *
 * Circle, line, points.
 */

class Point {
    double x;
    double y;

    public Point(double x, double y) {
        this.x = x;
        this.y = y;
    }

    @Override
    public String toString() {
        return "Point{" +
            "x=" + x +
            ", y=" + y +
            '}';
    }
}

class Circle {
    Point center;
    double radius;

    public Circle(Point center, double radius) {
        this.center = center;
        this.radius = radius;
    }

    @Override
    public String toString() {
        return "Circle{" +
            "center=" + center +
            ", radius=" + radius +
            '}';
    }
}

class Line {
    Point point;
    Point slope;

    public Line(Point point, Point slope) {
        this.point = point;
        this.slope = slope;
    }

    @Override
    public String toString() {
        return "Line{" +
            "point=" + point +
            ", slope=" + slope +
            '}';
    }
}

class GeoUtils {
    double error;
    double errorSqrt;
    double errorSq;
    double errorSix;

    /**
     * Constructs a geometry utility object. Point comparison is such that (a, b) and (a + e, b + e) are considered
     * the same. Other comparisons are scaled appropriately, e.g. to comparing two points on the unit circle.
     *
     * Hashing takes sqrt(error) by sqrt(error) sized blocks.
     * @param error The allowed error for comparisons
     */
    public GeoUtils(double error) {
        this.error = error;
        this.errorSqrt = Math.sqrt(error);
        this.errorSq = error * error;
        this.errorSix = errorSq * errorSq * errorSq;
    }

    // ========================================
    // Point operations
    // ========================================

    public double distanceSq(Point a, Point b) {
        double dx = a.x - b.x;
        double dy = a.y - b.y;
        return dx * dx + dy * dy;
    }

    public double distance(Point a, Point b) {
        return Math.sqrt(distanceSq(a, b));
    }

    public double lengthSq(Point a) {
        return a.x * a.x + a.y * a.y;
    }

    public Point sub(Point a, Point b) {
        return new Point(a.x - b.x, a.y - b.y);
    }

    public double dot(Point a, Point b) {
        return a.x * b.x + a.y * b.y;
    }

    // ========================================
    // Geometry operations
    // ========================================

    public Line line(Point a, Point b) {
        return new Line(a, sub(b, a));
    }

    public Circle circle(Point a, Point b) {
        return new Circle(a, distance(a, b));
    }

    public Point midpoint(Point a, Point b) {
        return new Point((a.x + b.x) / 2, (a.y + b.y) / 2);
    }


    // ========================================
    // Projection
    // ========================================

    public double projectionScalar(Point a, Point b) {
        return dot(a, b) / lengthSq(b);
    }

    public Point projection(Point a, Point b) {
        double scalar = projectionScalar(a, b);
        return new Point(scalar * b.x, scalar * b.y);
    }

    public Point rejection(Point a, Point b) {
        Point proj = projection(a, b);
        return sub(a, proj);
    }

    // ========================================
    // Contains
    // ========================================

    public boolean contains(Line l, Point p) {
        return parallel(sub(p, l.point), l.slope);
    }

    public boolean contains(Circle c, Point p) {
        return Math.abs(c.radius * c.radius - distanceSq(c.center, p)) < error * 2 * c.radius;
    }

    // ========================================
    // Utility
    // ========================================

    public long hash(Point p) {
        return ((long) (int) (p.x / errorSqrt)) * Integer.MAX_VALUE + (long) (int) (p.y / errorSqrt);
    }

    // ========================================
    // Equality
    // ========================================

    public boolean equals(double a, double b) {
        return Math.abs(a - b) < error;
    }

    public boolean equals(Point a, Point b) {
        return equals(a.x, b.x) && equals(a.y, b.y);
    }

    public boolean equals(Circle a, Circle b) {
        return equals(a.center, b.center) && equals(a.radius, b.radius);
    }

    public boolean equals(Line a, Line b) {
        return contains(a, b.point) && parallel(a.slope, b.slope);
    }

    public boolean parallel(Point a, Point b) {
        double cross = a.x * b.y - b.x * a.y;
        return cross * cross <= Math.max(errorSq * lengthSq(a) * lengthSq(b), errorSix);
    }

    public boolean parallel(Line a, Line b) {
        return parallel(a.slope, b.slope);
    }

    public boolean perpendicular(Line a, Line b) {
        double dot = dot(a.slope, b.slope);
        return dot * dot <= Math.max(errorSq * lengthSq(a.slope) * lengthSq(b.slope), errorSix);
    }

    // ========================================
    // Intersection
    // ========================================

    public Point[] intersect(Circle a, Circle b) {
        double d2 = distanceSq(a.center, b.center);
        double r0 = a.radius;
        double r1 = b.radius;
        if (d2 > (r0 + r1) * (r0 + r1) || d2 < (r0 - r1) * (r0 - r1)) {
            return new Point[0];
        }

        double baseScalar = (r0 * r0 - r1 * r1 + d2) / (2 * d2);

        double baseX = a.center.x + baseScalar * (b.center.x - a.center.x);
        double baseY = a.center.y + baseScalar * (b.center.y - a.center.y);
        double hScaledSq = r0 * r0 / d2 - baseScalar * baseScalar;
        double hScaled = hScaledSq > 0 ? Math.sqrt(hScaledSq) : 0;

        Point p1 = new Point(baseX + hScaled * (b.center.y - a.center.y), baseY - hScaled * (b.center.x - a.center.x));
        Point p2 = new Point(baseX - hScaled * (b.center.y - a.center.y), baseY + hScaled * (b.center.x - a.center.x));
        if (equals(p1, p2)) {
            return new Point[]{p1};
        }

        return new Point[]{p1, p2};
    }

    public Point[] intersect(Line a, Line b) {
        if (parallel(a, b)) {
            return new Point[0];
        }
        double denom = a.slope.y * b.slope.x - a.slope.x * b.slope.y;
        double scalar = (b.point.y - a.point.y) * b.slope.x - (b.point.x - a.point.x) * b.slope.y;
        return new Point[]{new Point(scalar / denom * a.slope.x + a.point.x, scalar / denom * a.slope.y + a.point.y)};
    }

    public Point[] intersect(Line a, Circle b) {
        return intersect(b, a);
    }

    public Point[] intersect(Circle a, Line b) {
        double dx = b.slope.x;
        double dy = b.slope.y;
        double bax = b.point.x - a.center.x;
        double bay = b.point.y - a.center.y;
        double A = dx * dx + dy * dy;
        double B = 2 * (dx * bax + dy * bay);
        double C = bax * bax + bay * bay - a.radius * a.radius;
        double det = B * B - 4 * A * C;
        if (det < 0) {
            return new Point[0];
        }

        double detSqrt = Math.sqrt(det);

        double t0 = (-B + detSqrt) / (2 * A);
        double t1 = (-B - detSqrt) / (2 * A);
        Point p1 = new Point(b.point.x + b.slope.x * t0, b.point.y + b.slope.y * t0);
        Point p2 = new Point(b.point.x + b.slope.x * t1, b.point.y + b.slope.y * t1);
        if (equals(p1, p2)) {
            return new Point[]{p1};
        }

        return new Point[]{p1, p2};
    }
}

class PointCount {
    Point p;
    int count;

    public PointCount(Point p, int count) {
        this.p = p;
        this.count = count;
    }
}

class State {
    private static class Operation {
        boolean circle;
        double a;
        double b;
        double c;
        double d;

        public Operation(Circle c) {
            this.circle = true;
            this.a = c.center.x;
            this.b = c.center.y;
            this.c = c.radius;
            this.d = c.radius;
        }

        public Operation(Line l) {
            this.circle = false;
            this.a = l.point.x;
            this.b = l.point.y;
            this.c = l.slope.x;
            this.d = l.slope.y;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Operation that = (Operation) o;

            if (circle != that.circle) return false;
            if (Double.compare(that.a, a) != 0) return false;
            if (Double.compare(that.b, b) != 0) return false;
            if (Double.compare(that.c, c) != 0) return false;
            return Double.compare(that.d, d) == 0;

        }

        @Override
        public int hashCode() {
            int result;
            long temp;
            result = (circle ? 1 : 0);
            temp = Double.doubleToLongBits(a);
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            temp = Double.doubleToLongBits(b);
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            temp = Double.doubleToLongBits(c);
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            temp = Double.doubleToLongBits(d);
            result = 31 * result + (int) (temp ^ (temp >>> 32));
            return result;
        }
    }

    GeoUtils geo;
    Circle[] circles;
    int cs;
    Line[] lines;
    int ls;
    Point[] points;
    int ps;

    HashMap<Long, PointCount> seen;

    // Same operations in a different order result in the same configuration
    long hash;

    // Operations on doubles are deterministic
    Random gen;
    HashMap<Operation, Long> zobrist;

    public State(GeoUtils geo, int maxDepth) {
        this.geo = geo;
        // TODO generalize later
        circles = new Circle[maxDepth + 10];
        lines = new Line[maxDepth + 10];
        points = new Point[maxDepth * maxDepth * 10];
        seen = new HashMap<>();
        cs = 0;
        ls = 0;
        ps = 0;
        hash = 0;
        gen = new Random(1337);
        zobrist = new HashMap<>();
    }

    public void add(Circle c) {
        pushCircle(c);
    }

    public void add(Line l) {
        pushLine(l);
    }

    public boolean add(Point p) {
        // Only check the hash box for the point.
        // There is a very small chance we add an already added point.
        long hash = geo.hash(p);
        if (seen.containsKey(hash)) {
            PointCount pointCount = seen.get(hash);
            if (pointCount.p != null && geo.equals(pointCount.p, p)) {
                return false;
            }

            // Do a full scan
            for (int i = 0; i < ps; i++) {
                if (geo.equals(points[i], p)) {
                    return false;
                }
            }

            if (pointCount.p == null) {
                pointCount.p = p;
            }
            pointCount.count++;
        } else {
            seen.put(hash, new PointCount(p, 1));
        }
        points[ps++] = p;
        return true;
    }

    public boolean containsPoint(Point p) {
        long hash = geo.hash(p);
        if (seen.containsKey(hash)) {
            PointCount pointCount = seen.get(hash);
            if (pointCount.p != null && geo.equals(pointCount.p, p)) {
                return true;
            }

            // Do a full scan
            for (int i = 0; i < ps; i++) {
                if (geo.equals(points[i], p)) {
                    return true;
                }
            }
        }

        return false;
    }

    public void removePoint(Point p) {
        long hash = geo.hash(p);
        PointCount pointCount = seen.get(hash);
        if (pointCount.p == p) {
            pointCount.p = null;
        }

        pointCount.count--;
        if (pointCount.count == 0) {
            seen.remove(hash);
        }
    }

    public boolean containsCircle(Circle circle) {
        for (int i = 0; i < cs; i++) {
            if (geo.equals(circle, circles[i])) {
                return true;
            }
        }

        return false;
    }

    public Circle getCircle(int a, int b) {
        return geo.circle(points[a], points[b]);
    }

    public int pushCircle(Circle circle) {
        int result = ps;

        for (int i = 0; i < cs; i++) {
            Point[] intersect = geo.intersect(circle, circles[i]);
            for (int j = 0; j < intersect.length; j++) {
                add(intersect[j]);
            }
        }

        for (int i = 0; i < ls; i++) {
            Point[] intersect = geo.intersect(circle, lines[i]);
            for (int j = 0; j < intersect.length; j++) {
                add(intersect[j]);
            }
        }

        circles[cs++] = circle;
        Operation o = new Operation(circle);
        zobrist.putIfAbsent(o, gen.nextLong());
        hash ^= zobrist.get(o);
        return result;
    }

    public void popCircle(int c) {
        cs--;
        hash ^= zobrist.get(new Operation(circles[cs]));

        for (int i = c; i < ps; i++) {
            removePoint(points[i]);
        }

        ps = c;
    }

    public Line getLine(int a, int b) {
        return geo.line(points[a], points[b]);
    }

    public boolean containsLine(Line line) {
        for (int i = 0; i < ls; i++) {
            if (geo.equals(line, lines[i])) {
                return true;
            }
        }

        return false;
    }

    public int pushLine(Line line) {
        int result = ps;

        for (int i = 0; i < cs; i++) {
            Point[] intersect = geo.intersect(line, circles[i]);
            for (int j = 0; j < intersect.length; j++) {
                add(intersect[j]);
            }
        }

        for (int i = 0; i < ls; i++) {
            Point[] intersect = geo.intersect(line, lines[i]);
            for (int j = 0; j < intersect.length; j++) {
                add(intersect[j]);
            }
        }

        lines[ls++] = line;
        Operation o = new Operation(line);
        zobrist.putIfAbsent(o, gen.nextLong());
        hash ^= zobrist.get(o);
        return result;
    }

    public void popLine(int c) {
        ls--;
        hash ^= zobrist.get(new Operation(lines[ls]));

        for (int i = c; i < ps; i++) {
            removePoint(points[i]);
        }

        ps = c;
    }
}

interface Problem {
    State getSetUp(int depth);
    boolean isSolved(State s);
}

class ProblemSolver {
    static class Operation {
        enum Type {
            CIRCLE(true), LINE(false), PERP_BISECT(false), PERP(false), ANGLE_BISECT(false), PARALLEL(false), CONG_CIRCLE(true);

            boolean circle;

            Type(boolean circle) {
                this.circle = circle;
            }
        }

        Type type;
        int a;
        int b;

        public Operation(Type type, int a, int b) {
            this.type = type;
            this.a = a;
            this.b = b;
        }

        @Override
        public String toString() {
            return "Operation{" +
                "type=" + type +
                ", a=" + a +
                ", b=" + b +
                '}';
        }
    }

    class Solution {
        State state;
        ArrayList<Operation> history;
        long counter;

        public Solution(State state, ArrayList<Operation> history, long counter) {
            this.state = state;
            this.history = history;
            this.counter = counter;
        }
    }

    // global search state
    HashSet<Long> seen;
    State s;
    int maxDepth;
    Problem p;
    long counter;

    ProblemSolver() {}

    public Solution solve(Problem p, int depth) {
        seen = new HashSet<>();
        s = p.getSetUp(depth);
        maxDepth = depth;
        counter = 0;
        this.p = p;
        Solution found = dfs(0);
        if (found == null) {
            return new Solution(null, null, counter);
        }

        Collections.reverse(found.history);
        return found;
    }

    public Solution dfs(int depth) {
        if (depth == maxDepth) {
            counter++;
            return p.isSolved(s) ? new Solution(s, new ArrayList<>(), counter) : null;
        }

        if (seen.contains(s.hash)) {
            return null;
        }

        seen.add(s.hash);
        for (int i = 0; i < s.ps; i++) {
            for (int j = 0; j < s.ps; j++) {
                if (j == i) {
                    continue;
                }

                Circle circle = s.getCircle(i, j);
                if (s.containsCircle(circle)) {
                    continue;
                }

                int c = s.pushCircle(circle);
                Solution found = dfs(depth + 1);

                if (found != null) {
                    found.history.add(new Operation(Operation.Type.CIRCLE, i, j));
                    return found;
                }
                s.popCircle(c);
            }
        }

        // Line constructions are symmetric
        for (int i = 0; i < s.ps; i++) {
            for (int j = i + 1; j < s.ps; j++) {
                Line line = s.getLine(i, j);
                if (s.containsLine(line)) {
                    continue;
                }

                int c = s.pushLine(line);
                Solution found = dfs(depth + 1);

                if (found != null) {
                    found.history.add(new Operation(Operation.Type.LINE, i, j));
                    return found;
                }
                s.popLine(c);
            }
        }

        return null;
    }
}

public class EuclideaSolver {
    public static boolean verifyCircleIntersect() {
        GeoUtils geo = new GeoUtils(0.0000001);
        Point p = new Point(0.4, 0.5);
        Point q = new Point(0.32, 0.65);
        Circle a = geo.circle(p, q);
        Circle b = geo.circle(q, p);
        Point[] intersect = geo.intersect(a, b);
        if (intersect.length != 2) {
            return false;
        }

        Point r = intersect[0];
        Circle c = geo.circle(r, q);
        Point[] intersectP = geo.intersect(c, b);
        if (intersectP.length != 2) {
            return false;
        }

        return geo.equals(intersectP[0], p) || geo.equals(intersectP[1], p);
    }

    public static boolean verifyLineIntersect() {
        GeoUtils geo = new GeoUtils(0.0000001);
        Point p = new Point(0.48, 0.95);
        Point q = new Point(0.34, 0.62);
        Circle a = geo.circle(p, q);
        Circle b = geo.circle(q, p);
        Point[] intersect = geo.intersect(a, b);
        if (intersect.length != 2) {
            return false;
        }

        Line l = geo.line(p, q);
        Line m = geo.line(intersect[0], intersect[1]);
        Point[] mid = geo.intersect(l, m);
        if (mid.length != 1) {
            return false;
        }

        return geo.equals(mid[0], geo.midpoint(p, q));
    }

    public static boolean verifyCircleLineIntersect() {
        GeoUtils geo = new GeoUtils(0.0000001);
        Point p = new Point(0.47, 0.856);
        Line l = new Line(p, new Point(0.42, 0.67));
        Point c = new Point(0.65, 0.537);
        Circle cp = geo.circle(c, p);

        Point[] intersect = geo.intersect(cp, l);
        if (intersect.length != 2) {
            return false;
        }

        Point q = intersect[0];
        if (geo.equals(p, q)) {
            q = intersect[1];
        }

        Line m = geo.line(q, c);
        Point[] top = geo.intersect(cp, m);
        if (top.length != 2) {
            return false;
        }

        Point t = top[0];
        if (geo.equals(t, p)) {
            t = top[1];
        }

        return geo.perpendicular(geo.line(p, t), l);
    }

    public static void verify() {
        System.out.println("VerifyCircleIntersect: " + verifyCircleIntersect());
        System.out.println("VerifyLineIntersect: " + verifyLineIntersect());
        System.out.println("VerifyCircleLineIntersect: " + verifyCircleLineIntersect());
    }

    public static void draw(BufferedImage image, Line l, Point origin, double scale, Color c) {
        Graphics graphics = image.getGraphics();
        graphics.setColor(c);

        int width = image.getWidth();
        int height = image.getHeight();
        int centerX = (int) ((l.point.x - origin.x) * scale + width / 2);
        int centerY = (int) ((l.point.y - origin.y) * scale + height / 2);
        int slopeX = (int)(l.slope.x * scale * 10);
        int slopeY = (int)(l.slope.y * scale * 10);
        graphics.drawLine(centerX - slopeX, centerY - slopeY, centerX + slopeX, centerY + slopeY);
    }

    public static void draw(BufferedImage image, Circle c, Point origin, double scale, Color d) {
        Graphics graphics = image.getGraphics();
        graphics.setColor(d);
        int width = image.getWidth();
        int height = image.getHeight();
        int centerX = (int) ((c.center.x - origin.x) * scale + width / 2);
        int centerY = (int) ((c.center.y - origin.y) * scale + height / 2);
        int radius = (int) (c.radius * scale);
        graphics.drawOval(centerX - radius, centerY - radius, radius * 2, radius * 2);
    }

    public static void showImage(BufferedImage image) {
        JFrame frame = new JFrame();
        frame.getContentPane().add(new JLabel(new ImageIcon(image)));
        frame.pack();
        frame.setVisible(true);
        frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
    }

    public static void drawAndShow(BufferedImage image, ProblemSolver.Solution solution) {
        State state = solution.state;
        int j = 0;
        int k = 0;

        Color[] colors = new Color[]{
            Color.RED, Color.ORANGE, Color.YELLOW, Color.GREEN, Color.CYAN, Color.BLUE, Color.MAGENTA,
        };

        int c = 0;
        int l = 0;
        for (ProblemSolver.Operation o : solution.history) {
            c += o.type.circle ? 1 : 0;
            l += o.type.circle ? 0 : 1;
        }

        j = state.cs - c;
        k = state.ls - l;
        for (int i = 0; i < j; i++) {
            draw(image, state.circles[i], new Point(0, 0), 100.0, Color.BLACK);
        }

        for (int i = 0; i < k; i++) {
            draw(image, state.lines[i], new Point(0, 0), 100.0, Color.BLACK);
        }

        for (int i = 0; i < solution.history.size(); i++) {
            ProblemSolver.Operation o = solution.history.get(i);

            Color color = i < colors.length ? colors[i] : Color.PINK;
            if (o.type.circle) {
                draw(image, state.circles[j++], new Point(0, 0), 100.0, color);
            } else {
                draw(image, state.lines[k++], new Point(0, 0), 100.0, color);
            }
        }

        showImage(image);
    }

    public static void main(String[] args) {
        verify();
        long start = System.currentTimeMillis();
        ProblemSolver solver = new ProblemSolver();
        ProblemSolver.Solution solution = solver.solve(new Problem9_1(), 5);
        System.out.println(System.currentTimeMillis() - start);
        System.out.println(solution.counter);

        if (solution.state == null) {
            return;
        }

        System.out.println(solution.history);
//        for (int i = 0; i < solution.state.ps; i++) {
//            System.out.println(solution.state.points[i]);
//        }
        BufferedImage image = new BufferedImage(600, 800, BufferedImage.TYPE_INT_ARGB);
        drawAndShow(image, solution);
    }
}

class Problem9_1 implements Problem {
    Point p = new Point(0, 0);
    Circle c = new Circle(new Point(0, 1), 0.32851489202);
    GeoUtils geo = new GeoUtils(1E-6);

    @Override
    public State getSetUp(int depth) {
        State state = new State(geo, depth);
        state.add(c);
        state.add(p);
        state.add(c.center);
        return state;
    }

    @Override
    public boolean isSolved(State s) {
        if (s.ls <= 0) {
            return false;
        }

        Line line = s.lines[s.ls - 1];
        if (!geo.contains(line, p)) {
            return false;
        }

        return geo.intersect(line, c).length == 1;
    }
}

class Problem9_9 implements Problem {
    Point p = new Point(0, 0);
    Point q = new Point(0, 1);
    Point goal = new Point(0, 1.0/6);
    GeoUtils geo = new GeoUtils(1E-6);
    Line l = geo.line(p, q);

    @Override
    public State getSetUp(int depth) {
        State state = new State(geo, depth);
        state.add(p);
        state.add(q);
        state.add(l);
        return state;
    }

    @Override
    public boolean isSolved(State s) {
        return s.containsPoint(goal);
    }
}
