#include <algorithm>
#include <array>
#include <cstdint>
#include <deque>
#include <iostream>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

// constexpr char input[] = "FlipSortMeZqjkYXWVUDCBghAN";
// constexpr char input[] = "FipSoUkYLDtMeZqjVhxCBgrwAN";
// constexpr char input[] = "okYLjDVtbpeZqFiuahxMgrwNSC";
constexpr char input[] = "abcdefgh";
// constexpr char input[] = "ABCDEF";
// constexpr char input[] = "ABCDEFGH"; // looks similar to pop count algorithm
// constexpr char input[] = "ABCDEFGHIJ"; // takes a long time
// constexpr char input[] = "AfBeCd";
// constexpr char input[] = "FlipSortMeZqjYBgAN";
// constexpr char input[] = "FlipSortMeZANJQG";
// constexpr char input[] = "FghqjklipnabcduvwxyzEmTSor";
// constexpr char input[] = "FgqjipnazEmTor";
// constexpr char input[] = "FghqjklipnabcduvwxyzEmTSRO";
// constexpr char input[] = "FEmnorstZYXWVUDCBghAKJILpq";
// constexpr char input[] = "FmrZDgAKLp"; // THIS ONE = 12 + 12 = 24.
// constexpr char input[] = "FgqlipnauEmTO";
// constexpr char input[] = "QFlipSortaMeZ";
// constexpr char input[] = "DfbecjavneAREJVLre";
// constexpr char input[] = "Sort";

constexpr size_t kMaxSize = sizeof(input) - 1;

class Node {
 private: 
  array<int8_t, kMaxSize> l;
  size_t s = kMaxSize;

 public: 
  // Node() {}

  Node(size_t s) : s(s) {}

  constexpr int8_t &operator[](size_t i) { return l[i]; }

  constexpr const int8_t &operator[](size_t i) const { return l[i]; }

  constexpr size_t size() const noexcept { return s; }

  typedef int8_t* iterator;
  typedef const int8_t* const_iterator;
  iterator begin() { return l.begin(); }
  const_iterator begin() const { return l.begin(); }
  iterator end() { return l.begin() + s; }
  const_iterator end() const { return l.begin() + s; }

  Node flip(int k) const {
    Node n{s};
    for (int i = 0; i < k; i++) {
      n.l[i] = l[i];
    }

    for (int i = k; i < s; i++) {
      n.l[s - 1 + k - i] = -l[i];
    }

    return n;
  }

  bool canSplit() {
    return s <= kMaxSize;
  }

  // adds an "a" to the beginning, and flips
  Node split() {
    Node n{s + 1};
    for (int i = 0; i < s; i++) {
      n.l[s - 1 - i] = -(l[i] + (2 * (l[i] > 0) - 1));
    }

    n.l[s] = -1;
    return n;
  }

  // splits the letter kth position and flips
  Node split(int k) {
    Node n{s + 1};
    for (int i = 0; i < k; i++) {
      // bithack: +1 if greater and positive, -1 if greater and negative, 0 if less.
      n.l[i] = l[i] + ((l[i] > l[k]) + (l[i] > -l[k]) - 1);
    }

    for (int i = k + 1; i < s; i++) {
      n.l[s + k - i] = -(l[i] + ((l[i] > l[k]) + (l[i] > -l[k]) - 1));
    }

    n.l[k] = l[k] - (l[k] < 0);
    n.l[s] = -(l[k] + (l[k] >= 0));
    return n;
  }

  Node canonical() const {
    Node n{s};
    for (int i = 0; i < s; i++) {
      n.l[i] = abs(l[i]);
    }

    sort(n.begin(), n.end());
    return n;
  }

  Node negate() const {
    Node n{s};
    for (int i = 0; i < s; i++) {
      n.l[i] = -l[i];
    }

    return n;
  }

  bool operator==(const Node &n) const {
    return s == n.s && equal(l.begin(), l.begin() + s, n.l.begin());
  }

  // Usable for kMaxSize <= 13
  size_t toInt() const {
    size_t h = 0;
    for (int i = 0; i < s; i++) {
      h = h * (kMaxSize * 2) + l[i];
    };
    return h;
  }
};

bool operator!=(const Node& l, const Node& r) {
  return !(l == r);
}

ostream &operator<<(ostream &os, const Node &n) {
  os << "[";
  for (int k : n) {
    os << k << ",";
  }
  os << "]";
  return os;
}

Node parse(const std::string &s) {
  Node n{s.size()};
  for (int i = 0; i < n.size(); i++) {
    if (s[i] >= 'a') {
      n[i] = s[i] - 'a' + 1;
    } else {
      n[i] = -(s[i] - 'A' + 1);
    }
  }
  return n;
}

string unparse(Node n) {
  string s;
  for (int i = 0; i < n.size(); i++) {
    if (n[i] >= 0) {
      s += 'a' + n[i] - 1;
    } else {
      s += 'A' - n[i] - 1;
    }
  }
  return s;
}

namespace std {
template<> struct hash<Node> {
  size_t operator()(const Node &n) const noexcept {
    size_t h = 0;
    for (int i = 0; i < n.size(); i++) {
      h = ((h << 5) + h) + n[i];
    };
    return h;
  }
};
} // namespace std

struct BlockMap {
  std::array<int8_t, 256> next = {0};
  BlockMap(Node n) {
    for (int i = 1; i < n.size(); i++) {
      next[n[i - 1] + 64] = n[i];
      next[-n[i] + 64] = -n[i - 1];
    }

    next[64] = n[0];
  }

  constexpr int8_t &operator[](int l) { return next[l + 64]; }

  constexpr const int8_t &operator[](int l) const { return next[l + 64]; }

  bool isJoined(const Node& u, int i) {
    return (i > 0 && (*this)[u[i - 1]] == u[i]) ||
           (i == 0 && (*this)[0] == u[i]);
  }

  bool isJoin(const Node& u, int i) {
    return (i > 0 && (*this)[u[i - 1]] == -u[u.size() - 1]) ||
           (i == 0 && (*this)[0] == -u[u.size() - 1]);
  }
};

bool isJoin(const Node &u, const Node &t, int i) {
  BlockMap next{t};
  return next.isJoined(u, i);
}

struct ScoredNode {
  Node n;
  int score;
  int id = maxId++;

  static int maxId;
  bool operator<(const ScoredNode &sn) const {
    return score > sn.score || (score == sn.score && id < sn.id);
  }
};

int ScoredNode::maxId = 0;

vector<int> astar(Node s, Node t) {
  if (s == t) {
    return {};
  }

  priority_queue<ScoredNode> queue;
  unordered_map<Node, int> prev;
  unordered_map<Node, int> dist;
  unordered_set<Node> closed;

  queue.push({s, 0});
  prev.emplace(s, -1);
  dist.emplace(s, 0);

  BlockMap nextT{t};
  BlockMap nextS{s};

  int lastDist = 0;
  int count = 0;
  while (!queue.empty()) {
    ScoredNode us = queue.top();
    Node u = us.n;

    if (count % 1024 == 0 && dist[u] > lastDist) {
      lastDist = dist[u];
      cout << "Finished " << lastDist << ", " << count << ", " << queue.size()
           << ", " << us.score << endl;
    }
    count++;

    if (u == t) {
      cout << "Finished " << dist[u] << ", " << count << ", " << queue.size()
           << ", " << us.score << endl;
      vector<int> flips;
      Node v = t;
      while (v != s) {
        int f = prev[v];
        flips.push_back(f);
        v = v.flip(f);
      }

      reverse(flips.begin(), flips.end());
      return flips;
    }
    queue.pop();

    if (closed.find(u) != closed.end()) {
      continue;
    }
    closed.emplace(u);

    for (int i = 0; i < u.size(); i++) {
      // Don't cut joined
      if (nextT.isJoined(u, i)) {
        continue;
      }

      // Don't join cut
      if (nextS.isJoin(u, i)) {
        continue;
      }

      Node v = u.flip(i);
      if (prev.find(v) != prev.end()) {
        continue;
      }

      int d = dist[u] + 1;
      if (dist.find(v) == dist.end() || d < dist[v]) {
        prev[v] = i;
        dist[v] = d;
        bool joined = nextT.isJoined(v, i);
        bool cut = nextS.isJoined(u, i);
        int score = us.score + 1 - joined;
        queue.push({v, score});
      }
    }
  }

  throw runtime_error("no path");
}

vector<int> bdbfs(Node s, Node t) {
  BlockMap smap{s};
  BlockMap tmap{t};

  deque<Node> sb;
  deque<Node> tb;
  sb.emplace_back(s);
  tb.emplace_back(t);

  if (s == t) {
    return {};
  }

  unordered_map<Node, int> sseen;
  sseen.emplace(s, -1);
  unordered_map<Node, int> tseen;
  tseen.emplace(t, -1);
  for (int i = 0;; i++) {
    {
      size_t size = sb.size();
      for (int j = 0; j < size; j++) {
        Node u = sb.front();
        sb.pop_front();
        for (int k = 0; k < u.size(); k++) {
          if (tmap.isJoined(u, k)) {
            continue;
          }

          Node v = u.flip(k);
          if (sseen.find(v) != sseen.end()) {
            continue;
          }

          if (tseen.find(v) != tseen.end()) {
            vector<int> flips;
            flips.resize(i * 2 + 1);
            Node p = v;
            int f = k;
            for (int a = i; a >= 0; a--) {
              flips[a] = f;
              p = p.flip(f);
              f = sseen[p];
            }
            p = v;
            f = tseen[p];
            for (int a = i + 1; a < 2 * i + 1; a++) {
              flips[a] = f;
              p = p.flip(f);
              f = tseen[p];
            }
            return flips;
          }

          sseen.emplace(v, k);
          sb.push_back(v);
        }
      }
    }
    {
      size_t size = tb.size();
      for (int j = 0; j < size; j++) {
        Node u = tb.front();
        tb.pop_front();

        for (int k = 0; k < u.size(); k++) {
          if (smap.isJoined(u, k)) {
            continue;
          }

          Node v = u.flip(k);
          if (tseen.find(v) != tseen.end()) {
            continue;
          }

          if (sseen.find(v) != sseen.end()) {
            vector<int> flips;
            flips.resize(i * 2 + 2);
            Node p = v;
            int f = k;
            for (int a = i + 1; a < 2 * i + 2; a++) {
              flips[a] = f;
              p = p.flip(f);
              f = tseen[p];
            }

            p = v;
            f = sseen[p];
            for (int a = i; a >= 0; a--) {
              flips[a] = f;
              p = p.flip(f);
              f = sseen[p];
            }
            return flips;
          }

          tseen.emplace(v, k);
          tb.push_back(v);
        }
      }
    }
    cout << "finished " << i << endl;
    cout << sb.size() << ", " << tb.size() << endl;
  }
}

// // return depth
int bfs(Node s, int limit) {
  // 1, 2, 3, 4, 5,  6,  7,  8,  9
  // 0, 2, 4, 6, 8, 10, 11, 13, 14
  BlockMap smap{s.canonical()};
  BlockMap clans{s.negate()};

  if (limit > kMaxSize) {
    throw runtime_error("limit > kMaxSize");
  }

  vector<deque<Node>> sbs;
  vector<size_t> sizes;
  sbs.resize(limit + 1);
  sizes.resize(limit + 1);
  sbs[0].emplace_back(Node{0});

  unordered_set<size_t> seen;
  seen.emplace(s.toInt());
  for (int i = 0;; i++) {
    cout << "Boundary " << i << ":\t";
    bool done = true;
    for (int m = 0; m < sbs.size(); m++) {
      sizes[m] = sbs[m].size();
      cout << sizes[m] << "\t";
      if (sizes[m] > 0) {
        done = false;
      }
    }
    cout << endl;

    for (int m = 0; m < sbs.size(); m++) {
      if (sbs[m].empty()) {
        continue;
      }
      // if (i < m * 3 / 2 + 2) {
        continue;
      // }

      for (auto& node : sbs[m]) {
        cout << unparse(node) << ", ";
      }
      cout << endl;
    }

    if (done) {
      return i - 1;
    }

    for (int m = 0; m <= limit; m++) {
      size_t size = sizes[m];
      auto& sb = sbs[m];
      for (int j = 0; j < size; j++) {
        Node u = sb.front();
        sb.pop_front();
        
        // flips
        for (int k = 0; k < u.size(); k++) {
          // // Only check things without clans
          // if (clans.isJoin(u, k)) {
          //   continue;
          // }
          if (smap.isJoin(u, k)) {
            continue;
          }

          Node v = u.flip(k);
          size_t vNum = v.toInt();
          if (seen.find(vNum) != seen.end()) {
            continue;
          }

          seen.emplace(vNum);
          sb.push_back(v);
        }

        // splits
        if (m == limit) {
          continue;
        }

        for (int k = 0; k < u.size(); k++) {
          Node v = u.split(k);
          size_t vNum = v.toInt();
          if (seen.find(vNum) != seen.end()) {
            continue;
          }

          seen.emplace(vNum);
          sbs[m + 1].push_back(v);
        }

        {
          Node v = u.split();
          size_t vNum = v.toInt();
          if (seen.find(vNum) != seen.end()) {
            continue;
          }

          seen.emplace(vNum);
          sbs[m + 1].push_back(v);
        }
      }
    }
  }
  return 0;
}

void mainShortest() {
  auto start = parse(input);
  auto fin = start.canonical();
  cout << unparse(start) << ", " << unparse(fin) << endl;
  auto path = astar(start, fin);
  // auto path = bdbfs(start, fin);

  cout << unparse(start) << endl;
  auto n = start;
  for (int i : path) {
    bool cut = isJoin(n, start, i);
    n = n.flip(i);
    bool join = isJoin(n, fin, i);
    cout << i << ", " << join << ", " << cut << endl;
    cout << unparse(n) << endl;
  }
  cout << path.size() << " steps" << endl;
}

void mainFill() {
  auto start = parse(input).canonical();
  cout << unparse(start) << endl;
  auto len = bfs(start, start.size());
  cout << len << " steps" << endl;
}

int main() {
  mainFill();
  // mainShortest();
  return 0;
}
