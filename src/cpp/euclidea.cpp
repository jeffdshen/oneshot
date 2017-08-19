#include <assert.h>
#include <math.h>
#include <stdint.h>

#include <algorithm>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <memory>
#include <random>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "SpookyV2.h"

using namespace std;

struct Point {
  double x;
  double y;

  bool operator==(const Point& p) const {
    return x == p.x && y == p.y;
  }
};

Point operator-(const Point& p, const Point& q) {
  return {p.x - q.x, p.y - q.y};
}

Point operator+(const Point& p, const Point& q) {
  return {p.x + q.x, p.y + q.y};
}

Point operator*(const Point& p, double scalar) {
  return {p.x * scalar, p.y * scalar};
}

Point operator*(double scalar, const Point& p) {
  return p * scalar;
}

Point operator/(const Point& p, double scalar) {
  return {p.x / scalar, p.y / scalar};
}

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

  double length(const Point& a) const {
    return sqrt(lengthSq(a));
  }

  Point rot90(const Point& a) const {
    return {-a.y, a.x};
  }

  double dot(const Point& a, const Point& b) const {
    return a.x * b.x + a.y * b.y;
  }

  // ========================================
  // Geometry operations
  // ========================================

  Line line(const Point& a, const Point& b) const {
    return {a, b - a};
  }

  Circle circle(const Point& a, const Point& b) const {
    return {a, distance(a, b)};
  }

  Circle congCircle(const Point& a, const Point& b, const Point& c) const {
    return {a, distance(b, c)};
  }

  Point midpoint(const Point& a, const Point& b) const {
    return (a + b) * 0.5;
  }

  Line perpBisect(const Point& a, const Point& b) const {
    return {midpoint(a, b), rot90(b - a)};
  }

  Line angleBisect(const Point& a, const Point& b, const Point& c) const {
    double bd = distance(a, b);
    double cd = distance(a, c);
    return line(a, (cd * b + bd * c) / (bd + cd));
  }

  Line perpendicular(const Point& p, const Line& line) {
    return {projection(p - line.point, line.slope) + line.point, rot90(line.slope)};
  }

  Line parallel(const Point& p, const Line& line) {
    return {p, line.slope};
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
    return a - proj;
  }

  // ========================================
  // Contains
  // ========================================

  bool contains(const Line& l, const Point& p) const {
    return parallel(p - l.point, l.slope);
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
 public:
  struct SetupSizes {
    size_t circles;
    size_t lines;
    size_t points;
  };

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
      uint64_t buf[5];
      buf[0] = k.circle;
      memcpy(buf + 1, &k.a, sizeof(k.a));
      memcpy(buf + 2, &k.b, sizeof(k.b));
      memcpy(buf + 3, &k.c, sizeof(k.c));
      memcpy(buf + 4, &k.d, sizeof(k.d));

      uint64_t hash1 = 0;
      uint64_t hash2 = 0;
      SpookyHash::Hash128((char *)buf, sizeof(buf), &hash1, &hash2);
      return hash1;
    }
  };

  GeoUtils geo_;
  vector<Circle> circles_;
  vector<Line> lines_;
  vector<Point> points_;
  SetupSizes setupSizes_ = {0, 0, 0};

  unordered_map<uint64_t, PointCount> seen_;

  // Same operations in a different order result in the same configuration
  uint64_t hash_;

  // Operations on doubles are deterministic
  std::default_random_engine gen_;
  std::uniform_int_distribution<uint64_t> dist_;

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
    hash_ ^= OperationHash{}(o);
    return result;
  }

  void removePoint(const Point& p) {
    uint64_t hash = geo_.hash(p);
    PointCount& pointCount = seen_.find(hash)->second;

    if (pointCount.count == 1) {
      seen_.erase(hash);
      return;
    }

    // exact match
    if (pointCount.p.first && pointCount.p.second == p) {
      pointCount.p.first = false;
    }

    pointCount.count--;
  }

 public:
  State(double error) : geo_(error), hash_(0), gen_(1337) {}
  State(const State&) = default;

  template <typename A, typename ... Args>
  State& push(A a, Args... args) {
    push(a);
    push(args...);
    return *this;
  };

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
    hash_ ^= OperationHash{}(o);
    circles_.pop_back();

    for (int i = c; i < points_.size(); i++) {
      removePoint(points_[i]);
    }

    points_.resize(c);
  }

  void popLine(int c) {
    Operation o(lines_.back());
    hash_ ^= OperationHash{}(o);
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

  void start() {
    setupSizes_ = {
        circles_.size(), lines_.size(), points_.size()
    };
  }

  const SetupSizes& getSetupSizes() const {
    return setupSizes_;
  }

  uint64_t getHash() const {
    return hash_;
  }

  const GeoUtils& getGeo() const {
    return geo_;
  }

  // ========================================
  // Utilities
  // ========================================

  // would do optional(pair<Point, Point>), std::optional when T.T
  bool canMake(const Line& line) const {
    // O(n), but still better than looping through all pairs.
    int hits = 0;
    for (auto& point : points_) {
      if (geo_.contains(line, point)) {
        hits++;
        if (hits >= 2) {
          return true;
        }
      }
    }

    return false;
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

// Used for caching hashes
template<typename Key>
class Cache {
 private:
  size_t capacity_;
  list<Key> list_;
  unordered_map<Key, typename list<Key>::iterator> map_;

public:
  Cache(size_t capacity) : capacity_(capacity) {}

  void insert(const Key& key) {
    list_.emplace_front(key);

    auto it = map_.find(key);
    if (it != map_.end()) {
      list_.erase(it->second);
    }

    map_[key] = list_.begin();

    if (map_.size() == capacity_ / 100) {
      cout << "almost there" << endl;
    }

    if (map_.size() > capacity_) {
      cout << "hit capacity" << endl;
      map_.erase(list_.back());
      list_.pop_back();
    }
  }

  bool contains(const Key& key) {
    return map_.find(key) != map_.end();
  }

  size_t capacity() {
    return capacity_;
  }

  size_t size() {
    return map_.size();
  }
};

struct Operation {
  enum class Type {
    CIRCLE, LINE, PERP_BISECT, PERP, ANGLE_BISECT, PARALLEL, CONG_CIRCLE
  };

  Type type;
  vector<uint32_t> points;

  Operation(Type t, const vector<uint32_t>& p) : type(t), points(p) {}

  static constexpr bool isSymmetric(Operation::Type type, int i, int j) {
    switch (type) {
      case Type::LINE:
      case Type::PERP_BISECT:
        return true;
      // focal point is first
      case Type::ANGLE_BISECT:
      case Type::CONG_CIRCLE:
        return i > 0 && j > 0;
      case Type::PERP:
      case Type::PARALLEL:
      case Type::CIRCLE:
        return false;
    }
  }

  static constexpr size_t numPoints(Operation::Type type) {
    switch (type) {
      case Type::CIRCLE:
      case Type::LINE:
      case Type::PERP_BISECT:
        return 2;
      case Type::PARALLEL:
      case Type::PERP:
      case Type::CONG_CIRCLE:
      case Type::ANGLE_BISECT:
        return 3;
    }
  }

  static constexpr bool isCircle(Operation::Type type) {
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

  static Type getType(const std::string& s) {
    if (s == "CIRCLE") {
      return Type::CIRCLE;
    } else if (s == "CONG_CIRCLE") {
      return Type::CONG_CIRCLE;
    } else if (s == "LINE") {
      return Type::LINE;
    } else if (s == "PERP_BISECT") {
      return Type::PERP_BISECT;
    } else if (s == "PERP") {
      return Type::PERP;
    } else if (s == "ANGLE_BISECT") {
      return Type::ANGLE_BISECT;
    } else if (s == "PARALLEL") {
      return Type::PARALLEL;
    } else {
      return Type::CIRCLE;
    }
  }

  bool isCircle() const {
    return isCircle(type);
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

struct CircleTool {
  Circle operator()(const State& s, size_t a, size_t b) {
    return s.getGeo().circle(s.getPoints()[a], s.getPoints()[b]);
  }
};

struct LineTool {
  Line operator()(const State& s, size_t a, size_t b) {
    return s.getGeo().line(s.getPoints()[a], s.getPoints()[b]);
  }
};

struct PerpBisectTool {
  Line operator()(const State& s, size_t a, size_t b) {
    return s.getGeo().perpBisect(s.getPoints()[a], s.getPoints()[b]);
  }
};

struct CongCircleTool {
  Circle operator()(const State& s, size_t a, size_t b, size_t c) {
    return s.getGeo().congCircle(s.getPoints()[a], s.getPoints()[b], s.getPoints()[c]);
  }
};

template <class Problem>
class ProblemSolverInstance {
 private:
  Problem p_;
  State s_;
  Cache<uint64_t> seen_;
//  unordered_set<uint64_t> seen_;
  vector<Operation> history_;
  const vector<set<Operation::Type>> guide_;
  const int maxDepth_;
  uint64_t counter_ = 0;
  uint64_t rejected_ = 0;

 public:
  ProblemSolverInstance(int maxDepth) : s_(p_.getSetup()), seen_(100000000), maxDepth_(maxDepth) {}

  ProblemSolverInstance(vector<set<Operation::Type>> guide) :
      s_(p_.getSetup()), seen_(100000000), guide_(guide), maxDepth_(guide.size()) {}

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

  template <class Tool>
  bool loopPointsTriple(Tool tool, Operation::Type type, int depth) {
    const auto& points = s_.getPoints();
    for (uint32_t i = 0; i < points.size(); i++) {
      uint32_t start = Operation::isSymmetric(type, 0, 1) ? i + 1 : 0;
      for (uint32_t j = start; j < points.size(); j++) {
        if (i == j) {
          continue;
        }

        start = Operation::isSymmetric(type, 0, 2) ? i + 1 : 0;
        start = max(start, Operation::isSymmetric(type, 1, 2) ? j + 1 : 0);
        for (uint32_t k = start; k < points.size(); k++) {
          if (k == j || k == i) {
            continue;
          }

          auto object = tool(s_, i, j, k);
          if (s_.contains(object)) {
            continue;
          }

          int c = s_.push(object);
          bool solved = dfs(depth + 1);

          if (solved) {
            history_.emplace_back(type, vector<uint32_t>{i, j, k});
            return true;
          }

          if (Operation::isCircle(type)) {
            s_.popCircle(c);
          } else {
            s_.popLine(c);
          }
        }
      }
      if (depth <= 2) {
        cout << depth << " : " << i + 1 << "/" << points.size() << endl;
      }
    }

    return false;
  }

  template <class Tool>
  bool loopPoints(Tool tool, Operation::Type type, int depth) {
    const auto& points = s_.getPoints();
    for (uint32_t i = 0; i < points.size(); i++) {
      int start = Operation::isSymmetric(type, 0, 1) ? i + 1 : 0;
      for (uint32_t j = start; j < points.size(); j++) {
        if (i == j) {
          continue;
        }

        auto object = tool(s_, i, j);
        if (s_.contains(object)) {
          continue;
        }

        int c = s_.push(object);
        bool solved = dfs(depth + 1);

        if (solved) {
          history_.emplace_back(type, vector<uint32_t>{i, j});
          return true;
        }

        if (Operation::isCircle(type)) {
          s_.popCircle(c);
        } else {
          s_.popLine(c);
        }
      }
      if (depth <= 2) {
        cout << depth << ":" << i + 1 << "/" << points.size() << endl;
      }
    }

    return false;
  }

  bool dfs(int depth) {
    if (!p_.isSolvable(depth, s_)) {
      counter_++;
      if (counter_ % 1000000 == 0) {
        cout << counter_ << endl;
        cout << seen_.size() << endl;
      }
      return false;
    }

    if (depth == maxDepth_) {
      counter_++;
      if (counter_ % 1000000 == 0) {
        cout << counter_ << endl;
        cout << seen_.size() << endl;
      }
      if (p_.isSolved(s_)) {
        return true;
      } else {
        return false;
      }
    }

    if (seen_.contains(s_.getHash())) {
      rejected_++;
      if (rejected_ % 1000000 == 0) {
        cout << "rejected: " << rejected_ << endl;
      }
//    if (seen_.find(s_.getHash()) != seen_.end()) {
      return false;
    }

    seen_.insert(s_.getHash());

    if (depth < guide_.size()) {
      auto& guide = guide_[depth];
      for (auto& type : guide) {
        bool solved = false;
        switch (type) {
          case Operation::Type::CIRCLE:
            CircleTool circle;
            solved = loopPoints(circle, Operation::Type::CIRCLE, depth);
            break;
          case Operation::Type::LINE:
            LineTool line;
            solved = loopPoints(line, Operation::Type::LINE, depth);
            break;
          case Operation::Type::PERP_BISECT:
            PerpBisectTool pb;
            solved = loopPoints(pb, Operation::Type::PERP_BISECT, depth);
            break;
          case Operation::Type::CONG_CIRCLE:
            CongCircleTool cc;
            solved = loopPointsTriple(cc, Operation::Type::CONG_CIRCLE, depth);
            break;
        }

        if (solved) {
          return true;
        }
      }

      return false;
    }


    bool solved = false;
    CircleTool circle;
    solved = loopPoints(circle, Operation::Type::CIRCLE, depth);
    if (solved) {
      return true;
    }

    LineTool line;
    solved = loopPoints(line, Operation::Type::LINE, depth);
    if (solved) {
      return true;
    }

//    PerpBisectTool pb;
//    solved = loopPoints(pb, Operation::Type::PERP_BISECT, depth);
//    if (solved) {
//      return true;
//    }

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
      << "state : " << s.state << ", "
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

  template<typename Problem>
  SolverResult solve(const vector<set<Operation::Type>>& guide) {
    ProblemSolverInstance<Problem> instance{guide};
    return instance.solve();
  }

  template<typename Problem>
  State replay(const vector<Operation>& history) {
    Problem p;
    State s = p.getSetup();
    for (Operation o : history) {
      switch (o.type) {
        case Operation::Type::CIRCLE:
          CircleTool circle;
          s.push(circle(s, o.points.at(0), o.points.at(1)));
          break;
        case Operation::Type::LINE:
          LineTool line;
          s.push(line(s, o.points.at(0), o.points.at(1)));
          break;
        case Operation::Type::PERP_BISECT:
          PerpBisectTool pb;
          s.push(pb(s, o.points.at(0), o.points.at(1)));
          break;
        case Operation::Type::CONG_CIRCLE:
          CongCircleTool cc;
          s.push(cc(s, o.points.at(0), o.points.at(1), o.points.at(2)));
          break;
      }
    }

    return s;
  }
};

vector<Operation> read(istream& in) {
  vector<Operation> history;
  while (!in.eof()) {
    std::string op;
    in >> op;
    // Only process operations, skip objects
    if (op == "Operation") {
      std::string typeName;
      in >> typeName;
      Operation::Type type = Operation::getType(typeName);
      vector<uint32_t> points;
      for (size_t i = 0; i < Operation::numPoints(type); i++) {
        uint32_t next;
        in >> next;
        points.emplace_back(next);
      }
      history.emplace_back(type, points);
    }
  }

  return history;
}

void print(ostream& out, const Circle& circle) {
  out << "Object" << " " << "Circle" << " "
      << circle.center.x << " " << circle.center.y << " "
      << circle.radius << "\n";
}

void print(ostream& out, const Line& line) {
  out << "Object" << " " << "Line" << " "
      << line.point.x << " " << line.point.y  << " "
      << line.slope.x << " " << line.slope.y << "\n";
}

void print(ostream& out, const Operation& operation) {
  out << "Operation" << " " << operation.getTypeName() << " ";
  if (!operation.points.empty()) {
    out << operation.points[0];
  }

  for (size_t i = 1; i < operation.points.size(); i++) {
    out << " " << operation.points[i];
  }

  out << "\n";
}

void print(ostream& out, const Solution& solution) {
  const auto& state = solution.state;
  size_t c = state.getSetupSizes().circles;
  size_t l = state.getSetupSizes().lines;

  const auto& circles = state.getCircles();
  for (size_t i = 0; i < c; i++) {
    print(out, circles[i]);
  }

  const auto& lines = state.getLines();
  for (size_t i = 0; i < l; i++) {
    print(out, lines[i]);
  }

  for (const auto& o : solution.history) {
    print(out, o);
    if (o.isCircle()) {
      print(out, circles[c++]);
    } else {
      print(out, lines[l++]);
    }
  }

  for (size_t i = c; i < circles.size(); i++) {
    print(out, circles[i]);
  }

  for (size_t i = l; i < lines.size(); i++) {
    print(out, lines[i]);
  }

  out << flush;
}

void writeSvg(ostream& out, const Point& point, size_t num) {
  out << "<text " << "x=\"" << point.x << "\" "
      << "y=\"" << point.y << "\" "
      << "font-size=\"0.5%\">"
      << num
      << "</text>\n";
}

void writeSvg(ostream& out, const Circle& circle, const string& color) {
  out << "<circle " << "cx=\"" << circle.center.x << "\" "
      << "cy=\"" << circle.center.y << "\" "
      << "r=\"" << circle.radius << "\" "
      << "stroke=\"" << color << "\" "
      << "stroke-width=\"0.05%\" fill-opacity=\"0.0\"/>\n";
}

void writeSvg(ostream& out, const Line& line, const string& color) {
  out << "<line " << "x1=\"" << line.point.x - 100 * line.slope.x << "\" "
      << "y1=\"" << line.point.y - 100 * line.slope.y << "\" "
      << "x2=\"" << line.point.x + 100 * line.slope.x << "\" "
      << "y2=\"" << line.point.y + 100 * line.slope.y << "\" "
      << "stroke=\"" << color << "\" "
      << "stroke-width=\"0.05%\" fill-opacity=\"0.0\"/>\n";
}

void writeSvgHeader(ostream& out, const Solution& solution) {
  double minX = 0;
  double minY = 0;
  double maxX = 0;
  double maxY = 0;
  for (const auto& point : solution.state.getPoints()) {
    minX = min(point.x, minX);
    minY = min(point.y, minY);
    maxX = max(point.x, maxX);
    maxY = max(point.y, maxY);
  }

  double midX = (minX + maxX) / 2;
  double radX = (maxX - minX) / 2;
  double midY = (minY + maxY) / 2;
  double radY = (maxY - minY) / 2;
  minX = midX - radX * 1.2; maxX = midX + radX * 1.2;
  minY = midY - radY * 1.2; maxY = midY + radY * 1.2;

  out << "<svg viewBox=\""
      << minX << " " << minY << " " << maxX - minX << " " << maxY - minY
      << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";
}

void writeSvg(ostream& out, const Solution& solution) {
  const auto& state = solution.state;
  const auto& circles = state.getCircles();
  const auto& lines = state.getLines();

  size_t c = state.getSetupSizes().circles;
  size_t l = state.getSetupSizes().lines;

  writeSvgHeader(out, solution);

  for (size_t i = 0; i < c; i++) {
    writeSvg(out, circles[i], "black");
  }

  for (size_t i = 0; i < l; i++) {
    writeSvg(out, lines[i], "black");
  }

  vector<string> colors = {
      "red", "hotpink", "orange", "green", "darkcyan", "blue", "violet", "purple", "brown"
  };

  size_t colorIndex = 0;
  for (const auto& o : solution.history) {
    string color = colorIndex < colors.size() ? colors[colorIndex] : "brown";
    colorIndex++;
    if (o.isCircle()) {
      writeSvg(out, circles[c++], color);
    } else {
      writeSvg(out, lines[l++], color);
    }
  }

  for (size_t i = c; i < circles.size(); i++) {
    writeSvg(out, circles[i], "silver");
  }

  for (size_t i = l; i < lines.size(); i++) {
    writeSvg(out, lines[i], "silver");
  }

  const auto& points = state.getPoints();
  for (size_t i = 0; i < points.size(); i++) {
    writeSvg(out, points[i], i);
  }
  out << "</svg>\n";
  out << flush;
}

// ========================================
// Problems
// ========================================

struct Problem {
  virtual bool isSolvable(int depth, State& s) {
    return true;
  }
};

struct Problem2_1 : Problem {
  GeoUtils geo{1E-7};
  Point p = {0, 0};
  Point a = {1, 0};
  Point b = {0.8482695305, -0.939458420};
  Line bisect = geo.angleBisect(p, a, b);

  State getSetup() {
    State state(1E-7);
    state.push(p, a, b);
    state.start();
    return state;
  }

  bool isSolved(const State& s) {
    return s.contains(bisect);
  }
};

struct Problem9_1 : Problem {
  Point p = {0, 0};
  Circle c = {{0, 1}, 0.32851489202};

  State getSetup() {
    State state(1E-7);
    state.push(c, p, c.center);
    state.start();
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

struct Problem9_9 : Problem {
  GeoUtils geo{1E-7};
  Point p{0, 0};
  Point q{0, 1};
  Point goal{0, 1.0 / 6};
  Line l = geo.line(p, q);

  State getSetup() {
    State state(1E-7);
    state.push(p, q, l);
    state.start();
    return state;
  }

  bool isSolved(const State& s) {
    return s.contains(goal);
  }
};

struct Problem9_2 : Problem {
  Circle c = {{0, 0}, 0.14293482043};
  Circle d = {{0, 1}, 0.32851489202};

  State getSetup() {
    State state(1E-7);
    state.push(c, d, c.center, d.center);
    state.start();
    return state;
  }

  bool isSolved(const State& s) {
    const auto& lines = s.getLines();
    if (lines.empty()) {
      return false;
    }

    const auto& geo = s.getGeo();
    const auto& line = lines.back();
    return geo.intersect(line, c).count == 1 && geo.intersect(line, d).count == 1;
  }
};

struct Problem10_1 : Problem {
  GeoUtils geo{1E-7};
  Point o{0, 0};
  Point a{0, 1};
  Point b{0, 1.254302934};
  Point c{0.432, 0.543};
  Point goal{b.y * c.x, b.y * c.y};
  Line l = geo.line(o, a);
  Line m = geo.line(o, c);

  State getSetup() {
    State state(1E-7);
    state.push(o, a, b, c, l, m);
    state.start();
    return state;
  }

  bool isSolved(const State& s) {
    return s.contains(goal);
  }
};

struct Problem10_2 : Problem {
  GeoUtils geo{1E-7};
  Point o{0, 0};
  Point a{0, 1};
  Point b{0, -0.34942294};
  Point c{0.432, sqrt(1 - 0.432 * 0.432)};
  Point goal{sqrt(-b.y) * c.x, sqrt(-b.y) * c.y};
  Line l = geo.line(o, a);
  Line m = geo.line(o, c);

  State getSetup() {
    State state(1E-7);
    state.push(o, a, b, l, m);
    state.start();
    return state;
  }

  bool isSolved(const State& s) {
    return s.contains(goal);
  }
};

struct Problem10_3 : Problem {
  GeoUtils geo{1E-7};
  Point o{0, 0};
  Point a{0, 1};
  Point goal{0, 0.61803398875};
  Line l = geo.line(o, a);

  State getSetup() {
    State state(1E-7);
    state.push(o, a, l);
    state.start();
    return state;
  }

  bool isSolved(const State& s) {
    return s.contains(goal);
  }
};

struct Problem10_5 : Problem {
  GeoUtils geo{1E-7};
  Point a{0, 0};
  Point b{1, 0};
  Point c{0, -1};
  Point d{1, -1};
  Point e{1.548186901, -1};
  Point p{0.652874589104, -1.43569032543};
  Line l = geo.line(a, b);
  Line m = geo.line(c, d);
  Line goal{p, {1, 0}};

  State getSetup() {
    State state(1E-7);
    state.push(l, m, c, d, e, p);
    state.start();
    return state;
  }

  bool isSolved(const State& s) {
    return s.getGeo().equals(s.getLines().back(), goal);
  }
};

struct Problem10_6LG : Problem {
  GeoUtils geo{1E-7};
  Point a{0, 0};
  Point b{1, 0};
  Point c{0.5426524491043, 0.3254264958193};
  Point p{0.7501967482255, 0.355652758493};
  Line bisect = geo.angleBisect(a, b, c);
  Line perp = geo.perpendicular(p, bisect);
  Line ab = geo.line(a, b);
  Line ac = geo.line(a, c);
  Circle goal = getGoal();

  State getSetup() {
    State state(1E-8);
    state.push(a, ab, ac, p);
    state.push(bisect); // step 1
    state.push(perp); // step 2
    state.start();
    return state;
  }

  Circle getGoal() {
    Line perp2 = geo.perpendicular(b, ab);
    Point o = geo.intersect(bisect, perp2).a;
    Circle circle = geo.circle(o, b);
    Line l = geo.line(a, p);
    Point p2 = geo.intersect(l, circle).a;
    Line para = geo.parallel(p, geo.line(o, p2));
    Point center = geo.intersect(para, bisect).a;
    return geo.circle(center, p);
  }

  bool isSolved(State& s) {
    if (s.getCircles().empty()) {
      return false;
    }

    for (auto& point : s.getPoints()) {
      auto& geo = s.getGeo();
      if (geo.equals(point, p)) {
        continue;
      }
      Line pb = geo.perpBisect(point, p);
      if (geo.equals(pb, bisect)) {
        continue;
      }

      if (geo.contains(pb, goal.center)) {
        s.push(pb); // step 5
        s.push(geo.circle(goal.center, p));  // step 6
        return true;
      }
    }

    return false;
  }
};

struct Problem11_3EG : Problem {
  GeoUtils geo{1E-7};
  Point a{0, 0};
  Point b{1, 0};
  Point c{0.5789492056831, -2.23584910493};
  Line al = geo.perpendicular(a, geo.angleBisect(a, b, c));
  Line bl = geo.perpendicular(b, geo.angleBisect(b, a, c));
  Line cl = geo.perpendicular(c, geo.angleBisect(c, a, b));

  State getSetup() {
    State state(1E-8);
    state.push(a, b, c);
    state.start();
    return state;
  }

  bool isSolved(State& state) {
    bool seen[] = {false, false, false};
    Line ll[] = {al, bl, cl};
    int pushed[] = {-1, -1, -1};
    size_t ps = 3;
    const auto& points = state.getPoints();
    const auto& geo = state.getGeo();

    for (size_t i = 0; i < ps; i++) {
      for (size_t j = 0; j < ps; j++) {
        if (seen[j]) {
          continue;
        }

        if (state.canMake(ll[j])) {
          seen[j] = true;
          pushed[i] = state.push(ll[j]); // step 8, 9, 10
          break;
        }
      }

      if (pushed[i] < 0) {
        for (int j = i - 1; j >= 0; j--) {
          state.popLine(pushed[j]);
        }

        return false;
      }
    }

    return true;
  }
};

struct Problem11_7EG : Problem {
  GeoUtils geo{1E-11};
  Point a{0, 0};
  Point b{1, 0};
  Line ab = geo.line(a, b);
  Point c{1.58433493, -1.44584910};
  Point d = c + Point{0.45849104958, 0.384939103948};
  Line cd = geo.line(c, d);
  Circle circle = geo.congCircle(b, c, d);
  Circle diam = geo.circle((a + b)/ 2, b);
  Point p = geo.intersect(circle, diam).a;
  Line l = geo.line(a, p);
  Line m = geo.line(b, p);

  State getSetup() {
    State state(1E-8);
    state.push(a, b, c, d);
    state.start();
    return state;
  }

  bool isSolvable(int depth, State& state) override {
    int d = depth;
    if (d == 6) {
      return state.canMake(l) || state.canMake(m);
    }

    if (d == 7) {
      const auto& geo = state.getGeo();
      const auto& lines = state.getLines();
      if (lines.empty()) {
        return false;
      }

      Line line = lines.back();
      if (!geo.equals(line, l) && !geo.equals(line, m)) {
        return false;
      }

      return true;
    }

    if (d == 8) {
      return state.canMake(l) && state.canMake(m);
    }

    return true;
  }

  bool isSolved(State& state) {
    int condition = 0;
//    GeoUtils geo = state.getGeo();
    for (auto& line : state.getLines()) {
      if (geo.equals(line, l)) {
        condition |= 1 << 0;
      } else if (geo.equals(line, m)) {
        condition |= 1 << 1;
      }

      if (condition >= (1 << 2) - 1) {
        return true;
      }
    }

    return false;
  }
};

// ========================================
// Main
// ========================================

using MainProblem = Problem11_7EG;

int main() {
//  Verifier v;
//  v.verify();
  ProblemSolver solver;

  auto ops = read(cin);
  State s = solver.replay<MainProblem>(ops);
  MainProblem p;
  bool solved = p.isSolved(s);
  cout << "Solved: " << solved << endl;
  Solution solution(s, ops);
  {
    print(cout, solution);
    ofstream file;
    file.open("solution.svg");
    writeSvg(file, solution);
    file.close();
    return 0;
  }

  {
    ofstream file;
    file.open("problem.svg");
    MainProblem problem;
    writeSvg(file, Solution(problem.getSetup(), {}));
    file.close();
  }

  auto start = chrono::steady_clock::now();
//  auto result = solver.solve<MainProblem>(5);
  auto result = solver.solve<MainProblem>(
      {
          {Operation::Type::LINE},
          {Operation::Type::CIRCLE},
          {Operation::Type::CIRCLE},
          {Operation::Type::CIRCLE},
          {Operation::Type::CIRCLE},
          {Operation::Type::CIRCLE},
          {Operation::Type::LINE},
          {Operation::Type::CIRCLE},
          {Operation::Type::LINE},
      }
  );
  auto end = chrono::steady_clock::now();
  cout << "Run took " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms" << endl;
  cout << "Checked " << result.counter << " end states" << endl;
  if (result.solution) {
    print(cout, *result.solution);
    ofstream file;
    file.open("euclidea.svg");
    writeSvg(file, *result.solution);
    file.close();
  }
}
