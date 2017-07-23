#include <assert.h>
#include <math.h>
#include <stdint.h>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <iterator>
#include <memory>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

struct Point {
  double x;
  double y;

  bool operator==(const Point& p) const {
    return x == p.x && y == p.y;
  }
};

std::ostream& operator<< (std::ostream& out, const Point& p) {
  out << "Point{ "
      << "x : " << p.x << ", "
      << "y : " << p.y << " }";
  return out ;
}

struct Circle {
  Point center;
  double radius;
};

std::ostream& operator<< (std::ostream& out, const Circle& c) {
  out << "Circle{ "
      << "center : " << c.center << ", "
      << "radius : " << c.radius << " }";
  return out ;
}

struct Line {
  Point point;
  Point slope;
};

std::ostream& operator<< (std::ostream& out, const Line& l) {
  out << "Line{ "
      << "point : " << l.point << ", "
      << "slope : " << l.slope << " }";
  return out ;
}

struct Intersection {
  int count;
  Point a;
  Point b;
};

class GeoUtils {
 private:
  const double error_;
  const double errorSqrt_;
  const double errorSq_;
  const double errorSix_;
 public:
  /**
  * Constructs a geometry utility object. Point comparison is such that (a, b) and (a + e, b + e) are considered
  * the same. Other comparisons are scaled appropriately, e.g. to comparing two points on the unit circle.
  *
  * Hashing takes sqrt(error) by sqrt(error) sized blocks.
  * @param error The allowed error for comparisons
  */
  GeoUtils(double error) :
      error_(error), errorSqrt_(sqrt(error)), errorSq_(error * error), errorSix_(errorSq_ * errorSq_ * errorSq_) {}

  GeoUtils(const GeoUtils&) = default;

  // ========================================
  // Point operations
  // ========================================

  double distanceSq(const Point& a, const Point& b) const {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return dx * dx + dy * dy;
  }

  double distance(const Point& a, const Point& b) const {
    return sqrt(distanceSq(a, b));
  }

  double lengthSq(const Point& a) const {
    return a.x * a.x + a.y * a.y;
  }

  Point sub(const Point& a, const Point& b) const {
    return {a.x - b.x, a.y - b.y};
  }

  double dot(const Point& a, const Point& b) const {
    return a.x * b.x + a.y * b.y;
  }

  // ========================================
  // Geometry operations
  // ========================================

  Line line(const Point& a, const Point& b) const {
    return {a, sub(b, a)};
  }

  Circle circle(const Point& a, const Point& b) const {
    return {a, distance(a, b)};
  }

  Point midpoint(const Point& a, const Point& b) const {
    return {(a.x + b.x) / 2, (a.y + b.y) / 2};
  }

  // ========================================
  // Projection
  // ========================================

  double projectionScalar(const Point& a, const Point& b) const {
    return dot(a, b) / lengthSq(b);
  }

  Point projection(const Point& a, const Point& b) const {
    double scalar = projectionScalar(a, b);
    return {scalar * b.x, scalar * b.y};
  }

  Point rejection(const Point& a, const Point& b) const {
    Point proj = projection(a, b);
    return sub(a, proj);
  }

  // ========================================
  // Contains
  // ========================================

  bool contains(const Line& l, const Point& p) const {
    return parallel(sub(p, l.point), l.slope);
  }

  bool contains(const Circle& c, const Point& p) const {
    return abs(c.radius * c.radius - distanceSq(c.center, p)) < error_ * 2 * c.radius;
  }

  // ========================================
  // Utility
  // ========================================

  uint64_t hash(const Point& p) const {
    return (static_cast<uint64_t>(static_cast<int64_t>(p.x / errorSqrt_)) & 0xFFFFFFFF) << 32 |
        static_cast<uint64_t>(static_cast<int64_t>(p.y / errorSqrt_)) & 0xFFFFFFFF;
  }

  // ========================================
  // Equality
  // ========================================

  bool equals(double a, double b) const {
    return abs(a - b) < error_;
  }

  bool equals(const Point& a, const Point& b) const {
    return equals(a.x, b.x) && equals(a.y, b.y);
  }

  bool equals(const Circle& a, const Circle& b) const {
    return equals(a.center, b.center) && equals(a.radius, b.radius);
  }

  bool equals(const Line& a, const Line& b) const {
    return contains(a, b.point) && parallel(a.slope, b.slope);
  }

  bool parallel(const Point& a, const Point& b) const {
    double cross = a.x * b.y - b.x * a.y;
    return cross * cross <= max(errorSq_ * lengthSq(a) * lengthSq(b), errorSix_);
  }

  bool parallel(const Line& a, const Line& b) const {
    return parallel(a.slope, b.slope);
  }

  bool perpendicular(const Line& a, const Line& b) const {
    double d = dot(a.slope, b.slope);
    return d * d <= max(errorSq_ * lengthSq(a.slope) * lengthSq(b.slope), errorSix_);
  }

  // ========================================
  // Intersection
  // ========================================

  Intersection intersect(const Circle& a, const Circle& b) const {
    double d2 = distanceSq(a.center, b.center);
    double r0 = a.radius;
    double r1 = b.radius;
    if (d2 > (r0 + r1) * (r0 + r1) || d2 < (r0 - r1) * (r0 - r1)) {
      return {0, {}, {}};
    }

    double baseScalar = (r0 * r0 - r1 * r1 + d2) / (2 * d2);
    double baseX = a.center.x + baseScalar * (b.center.x - a.center.x);
    double baseY = a.center.y + baseScalar * (b.center.y - a.center.y);
    double hScaledSq = r0 * r0 / d2 - baseScalar * baseScalar;
    double hScaled = hScaledSq > 0 ? sqrt(hScaledSq) : 0;

    Point p1{baseX + hScaled * (b.center.y - a.center.y), baseY - hScaled * (b.center.x - a.center.x)};
    Point p2{baseX - hScaled * (b.center.y - a.center.y), baseY + hScaled * (b.center.x - a.center.x)};
    if (equals(p1, p2)) {
      return {1, p1, {}};
    }

    return {2, p1, p2};
  }

  Intersection intersect(const Line& a, const Line& b) const {
    if (parallel(a, b)) {
      return {0, {}, {}};
    }
    double denom = a.slope.y * b.slope.x - a.slope.x * b.slope.y;
    double scalar = (b.point.y - a.point.y) * b.slope.x - (b.point.x - a.point.x) * b.slope.y;
    return {1, {scalar / denom * a.slope.x + a.point.x, scalar / denom * a.slope.y + a.point.y}, {}};
  }

  Intersection intersect(const Line& a, const Circle& b) const {
    return intersect(b, a);
  }

  Intersection intersect(const Circle& a, const Line& b) const {
    double dx = b.slope.x;
    double dy = b.slope.y;
    double bax = b.point.x - a.center.x;
    double bay = b.point.y - a.center.y;
    double A = dx * dx + dy * dy;
    double B = 2 * (dx * bax + dy * bay);
    double C = bax * bax + bay * bay - a.radius * a.radius;
    double det = B * B - 4 * A * C;
    if (det < 0) {
      return {0, {}, {}};
    }

    double detSqrt = sqrt(det);

    double t0 = (-B + detSqrt) / (2 * A);
    double t1 = (-B - detSqrt) / (2 * A);
    Point p1{b.point.x + b.slope.x * t0, b.point.y + b.slope.y * t0};
    Point p2{b.point.x + b.slope.x * t1, b.point.y + b.slope.y * t1};
    if (equals(p1, p2)) {
      return {1, p1, {}};
    }

    return {2, p1, p2};
  }
};

struct PointCount {
  pair<bool, Point> p; // "optional" field
  int count;
};

class State {
 private:
  struct Operation {
    bool circle;
    double a;
    double b;
    double c;
    double d;

    Operation(const Circle& o) :
        circle(true), a(o.center.x), b(o.center.y), c(o.radius), d(o.radius) {}

    Operation(const Line& l) : circle(false), a(l.point.x), b(l.point.y), c(l.slope.x), d(l.slope.y) {}

    bool operator==(const Operation& o) const {
      return circle == o.circle && a == o.a && b == o.b && c == o.c && d == o.d;
    }
  };

  struct OperationHash {
    size_t operator()(const Operation& k) const {
      std::hash<double> hasher;
      size_t result = k.circle;
      result = (result << 1) ^ hasher(k.a);
      result = (result << 1) ^ hasher(k.b);
      result = (result << 1) ^ hasher(k.c);
      result = (result << 1) ^ hasher(k.d);
      return result;
    }
  };

  GeoUtils geo_;
  vector<Circle> circles_;
  vector<Line> lines_;
  vector<Point> points_;

  unordered_map<uint64_t, PointCount> seen_;

  // Same operations in a different order result in the same configuration
  uint64_t hash_;

  // Operations on doubles are deterministic
  std::default_random_engine gen_;
  std::uniform_int_distribution<uint64_t> dist_;
  unordered_map<Operation, uint64_t, OperationHash> zobrist_;

 private:
  template<class A, class B>
  void addIntersects(const A& a, const vector<B>& bs) {
    for (auto& b : bs) {
      auto intersect = geo_.intersect(a, b);
      switch (intersect.count) {
        case 2:
          push(intersect.b);
          // fallthrough
        case 1:
          push(intersect.a);
      }
    }
  }

  template<class A>
  int addInternal(const A& a, vector<A>& as) {
    int result = points_.size();

    addIntersects(a, circles_);
    addIntersects(a, lines_);
    Operation o(a);

    as.emplace_back(a);
    zobrist_.emplace(o, dist_(gen_));
    hash_ ^= zobrist_.at(o);
    return result;
  }

  void removePoint(const Point& p) {
    uint64_t hash = geo_.hash(p);
    PointCount& pointCount = seen_.find(hash)->second;

    // exact match
    if (pointCount.p.first && pointCount.p.second == p) {
      pointCount.p.first = false;
    }

    pointCount.count--;
  }

 public:
  State(double error) : geo_(error), hash_(0), gen_(1337) {
    zobrist_.reserve(10000000);
  }
  State(const State&) = default;

  int push(const Circle& c) {
    return addInternal(c, circles_);
  }

  int push(const Line& l) {
    return addInternal(l, lines_);
  }

  bool push(const Point& p) {
    // Only check the hash box for the point.
    // There is a very small chance we add an already added point.
    uint64_t hash = geo_.hash(p);

    auto it = seen_.find(hash);
    if (it != seen_.end()) {
      PointCount& pointCount = it->second;
      if (pointCount.p.first && geo_.equals(pointCount.p.second, p)) {
        return false;
      }

      // Do a full scan
      for (auto& point : points_) {
        if (geo_.equals(point, p)) {
          return false;
        }
      }

      if (!pointCount.p.first) {
        pointCount.p.second = p;
      }
      pointCount.count++;
    } else {
      seen_.emplace(hash, PointCount{{true, p}, 1});
    }

    points_.emplace_back(p);
    return true;
  }

  bool contains(const Point& p) const {
    uint64_t hash = geo_.hash(p);
    auto it = seen_.find(hash);
    if (it != seen_.end()) {
      const PointCount& pointCount = it->second;
      if (pointCount.p.first && geo_.equals(pointCount.p.second, p)) {
        return true;
      }

      // Do a full scan
      for (auto& point : points_) {
        if (geo_.equals(point, p)) {
          return true;
        }
      }
    }

    return false;
  }

  bool contains(const Circle& c) const {
    for (auto& circle : circles_) {
      if (geo_.equals(circle, c)) {
        return true;
      }
    }

    return false;
  }

  bool contains(const Line& l) const {
    for (auto& line : lines_) {
      if (geo_.equals(line, l)) {
        return true;
      }
    }

    return false;
  }

  void popCircle(int c) {
    Operation o(circles_.back());
    hash_ ^= zobrist_.at(o);
    circles_.pop_back();

    for (int i = c; i < points_.size(); i++) {
      removePoint(points_[i]);
    }

    points_.resize(c);
  }

  void popLine(int c) {
    Operation o(lines_.back());
    hash_ ^= zobrist_.at(o);
    lines_.pop_back();

    for (int i = c; i < points_.size(); i++) {
      removePoint(points_[i]);
    }

    points_.resize(c);
  }

  const vector<Point>& getPoints() const {
    return points_;
  }

  const vector<Circle>& getCircles() const {
    return circles_;
  }

  const vector<Line>& getLines() const {
    return lines_;
  }

  uint64_t getHash() const {
    return hash_;
  }

  const GeoUtils& getGeo() const {
    return geo_;
  }
};

class Verifier {
 public:
  bool verifyCircleIntersect() {
    GeoUtils geo(1e-7);
    Point p{0.4, 0.5};
    Point q{0.32, 0.65};
    Circle a = geo.circle(p, q);
    Circle b = geo.circle(q, p);
    Intersection intersect = geo.intersect(a, b);
    if (intersect.count != 2) {
      return false;
    }

    Point r = intersect.a;
    Circle c = geo.circle(r, q);
    Intersection intersectP = geo.intersect(c, b);
    if (intersectP.count != 2) {
      return false;
    }

    return geo.equals(intersectP.a, p) || geo.equals(intersectP.b, p);
  }

  bool verifyLineIntersect() {
    GeoUtils geo(1e-7);
    Point p{0.48, 0.95};
    Point q{0.34, 0.62};
    Circle a = geo.circle(p, q);
    Circle b = geo.circle(q, p);
    Intersection intersect = geo.intersect(a, b);
    if (intersect.count != 2) {
      return false;
    }

    Line l = geo.line(p, q);
    Line m = geo.line(intersect.a, intersect.b);
    Intersection mid = geo.intersect(l, m);
    if (mid.count != 1) {
      return false;
    }

    return geo.equals(mid.a, geo.midpoint(p, q));
  }

  bool verifyCircleLineIntersect() {
    GeoUtils geo(1e-7);
    Point p{0.47, 0.856};
    Line l{p, {0.42, 0.67}};
    Point c {0.65, 0.537};
    Circle cp = geo.circle(c, p);

    Intersection intersect = geo.intersect(cp, l);
    if (intersect.count != 2) {
      return false;
    }

    Point q = intersect.a;
    if (geo.equals(p, q)) {
      q = intersect.b;
    }

    Line m = geo.line(q, c);
    Intersection top = geo.intersect(cp, m);
    if (top.count != 2) {
      return false;
    }

    Point t = top.a;
    if (geo.equals(t, p)) {
      t = top.b;
    }

    return geo.perpendicular(geo.line(p, t), l);
  }

 void verify() {
    cout << "VerifyCircleIntersect: " << verifyCircleIntersect() << endl;
    cout << "VerifyLineIntersect: " << verifyLineIntersect() << endl;
    cout << "VerifyCircleLineIntersect: " << verifyCircleLineIntersect() << endl;
  }
};

struct Operation {
  enum class Type {
    CIRCLE, LINE, PERP_BISECT, PERP, ANGLE_BISECT, PARALLEL, CONG_CIRCLE
  };

  Type type;
  vector<uint32_t> points;

  Operation(Type t, const vector<uint32_t>& p) : type(t), points(p) {}

  bool isCircle() {
    switch (type) {
      case Type::CIRCLE:
      case Type::CONG_CIRCLE:
        return true;
      case Type::LINE:
      case Type::PERP_BISECT:
      case Type::PERP:
      case Type::ANGLE_BISECT:
      case Type::PARALLEL:
        return false;
    }
  }

  std::string getTypeName() const {
    switch (type) {
      case Type::CIRCLE:
        return "CIRCLE";
      case Type::CONG_CIRCLE:
        return "CONG_CIRCLE";
      case Type::LINE:
        return "LINE";
      case Type::PERP_BISECT:
        return "PERP_BISECT";
      case Type::PERP:
        return "PERP";
      case Type::ANGLE_BISECT:
        return "ANGLE_BISECT";
      case Type::PARALLEL:
        return "PARALLEL";
    }
  }
};

struct Solution {
  State state;
  vector<Operation> history;

  Solution(const State& s, const vector<Operation>& h) : state(s), history(h) {}
};

struct SolverResult {
  unique_ptr<Solution> solution;
  long counter;
};

template <class Problem>
class ProblemSolverInstance {
 private:
  Problem p_;
  State s_;
  unordered_set<uint64_t> seen_;
  vector<Operation> history_;
  const int maxDepth_;
  uint64_t counter_;

 public:
  ProblemSolverInstance(int maxDepth) : s_(p_.getSetup()), maxDepth_(maxDepth), counter_(0) {}

  SolverResult solve() {
    SolverResult result;
    bool solved = dfs(0);

    if (solved) {
      result.solution = make_unique<Solution>(s_, history_);
      reverse(result.solution->history.begin(), result.solution->history.end());
    }

    result.counter = counter_;
    return result;
  }

  bool dfs(int depth) {
    if (depth == maxDepth_) {
      counter_++;
      if (counter_ % 1000000 == 0) {
        cout << counter_ << endl;
      }
      if (p_.isSolved(s_)) {
        return true;
      } else {
        return false;
      }
    }

    if (seen_.find(s_.getHash()) != seen_.end()) {
      return false;
    }

    seen_.insert(s_.getHash());

    const auto& points = s_.getPoints();
    for (uint32_t i = 0; i < points.size(); i++) {
      for (uint32_t j = 0; j < points.size(); j++) {
        if (i == j) {
          continue;
        }

        Circle circle = s_.getGeo().circle(points[i], points[j]);
        if (s_.contains(circle)) {
          continue;
        }

        int c = s_.push(circle);
        bool solved = dfs(depth + 1);

        if (solved) {
          history_.emplace_back(Operation::Type::CIRCLE, vector<uint32_t>{i, j});
          return true;
        }

        s_.popCircle(c);
      }
    }

    // Line constructions are symmetric
    for (uint32_t i = 0; i < points.size(); i++) {
      for (uint32_t j = i + 1; j < points.size(); j++) {
        Line line = s_.getGeo().line(points[i], points[j]);
        if (s_.contains(line)) {
          continue;
        }

        int c = s_.push(line);
        bool solved = dfs(depth + 1);

        if (solved) {
          history_.emplace_back(Operation::Type::LINE, vector<uint32_t>{i, j});
          return true;
        }

        s_.popLine(c);
      }
    }

    return false;
  }
};

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
  if ( !v.empty() ) {
    out << '[';
    copy(v.begin(), v.end(), ostream_iterator<T>(out, ", "));
    out << "\b\b]";
  }
  return out;
}


std::ostream& operator<< (std::ostream& out, const State& s) {
  out << "State{"
      << "circles : " << s.getCircles() << ", "
      << "lines : " << s.getLines() << ", "
      << "points : " << s.getPoints() << " }";
  return out;
}


std::ostream& operator<< (std::ostream& out, const Operation& o) {
  out << "Operation{"
      << "type : " << o.getTypeName() << ", "
      << "points : " << o.points << "}";
  return out;
}

std::ostream& operator<<(std::ostream& out, const Solution& s) {
  out << "Solution{"
      << "state : " << s.state << " }"
      << "history : " << s.history << " }";
  return out;
}

class ProblemSolver {
 public:
  template<typename Problem>
  SolverResult solve(int depth) {
    ProblemSolverInstance<Problem> instance{depth};
    return instance.solve();
  }
};

struct Problem9_1 {
  Point p = {0, 0};
  Circle c = {{0, 1}, 0.32851489202};

  State getSetup() {
    State state(1E-6);
    state.push(c);
    state.push(p);
    state.push(c.center);
    return state;
  }

  bool isSolved(const State& s) {
    const auto& lines = s.getLines();
    if (lines.empty()) {
      return false;
    }

    const auto& geo = s.getGeo();
    const auto& line = lines.back();
    if (!geo.contains(line, p)) {
      return false;
    }

    return geo.intersect(line, c).count == 1;
  }
};

struct Problem9_9 {
  GeoUtils geo{1E-6};
  Point p{0, 0};
  Point q{0, 1};
  Point goal{0, 1.0 / 6};
  Line l = geo.line(p, q);

  State getSetup() {
    State state(1E-6);
    state.push(p);
    state.push(q);
    state.push(l);
    return state;
  }

  bool isSolved(const State& s) {
    return s.contains(goal);
  }
};

int main() {
  Verifier v;
  v.verify();
  ProblemSolver solver;
  auto start = chrono::steady_clock::now();
  auto result = solver.solve<Problem9_1>(5);
  auto end = chrono::steady_clock::now();
  cout << chrono::duration_cast<chrono::milliseconds>(end - start).count() << endl;
  cout << result.counter << endl;
  if (result.solution) {
    cout << *(result.solution) << endl;
  }
}
