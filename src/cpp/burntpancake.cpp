#include <algorithm>
#include <array>
#include <cstdint>
#include <deque>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>

using namespace std;

constexpr char input[] = "FlipSortMeZqjkYXWVUDCBghAN";
// constexpr char input[] = "FipSoUkYLDtMeZqjVhxCBgrwAN";
// constexpr char input[] = "okYLjDVtbpeZqFiuahxMgrwNSC";
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

struct Node {
  array<int8_t, sizeof(input) - 1> l;

  constexpr int8_t& operator[](size_t i) {
    return l[i];
  }

  constexpr const int8_t& operator[](size_t i) const {
    return l[i];
  }

  constexpr size_t size() const noexcept {
    return l.size();
  }

  Node flip(int k) const {
    Node n;
    for (int i = 0; i < k ; i++) {
      n.l[i] = l[i];
    }

    for (int i = k; i < l.size(); i++) {
      n.l[i] = -l[l.size() - 1 + k - i];
    }

    return n;
  }

  Node canonical() const {
    Node n;
    for (int i = 0; i < l.size(); i++) {
      n.l[i] = abs(l[i]);
    }

    sort(n.l.begin(), n.l.end());
    return n;
  }

  bool operator==(const Node &n) const { return l == n.l; }

  bool operator!=(const Node &n) const { return l != n.l; }
};

ostream &operator<<(ostream &os, const Node &n) {
  os << "[";
  for (int k : n.l) {
    os << k << ",";
  }
  os << "]";
  return os;
}

Node parse(const std::string &s) {
  Node n;
  for (int i = 0; i < n.l.size(); i++) {
    if (s[i] >= 'a') {
      n.l[i] = s[i] - 'a' + 1;
    } else {
      n.l[i] = -(s[i] - 'A' + 1);
    }
  }
  return n;
}

string unparse(Node n) {
  string s;
  for (int i = 0; i < n.l.size(); i++) {
    if (n.l[i] >= 0) {
      s += 'a' + n.l[i] - 1;
    } else {
      s += 'A' - n.l[i] - 1;
    }
  }
  return s;
}

namespace std {
template <> struct hash<Node> {
  size_t operator()(const Node &s) const noexcept {
    size_t h = 0;
    for (int i = 0; i < s.l.size(); i++) {
      h  = ((h << 5) + h) + s.l[i];
    };
    return h;
  }
};
} // namespace std

struct NextMap {
  std::array<int8_t, 256> next = {0};
  NextMap(Node n) {
    for (int i = 1; i < n.l.size(); i++) {
      next[n.l[i-1] + 64] = n.l[i];
      next[-n.l[i] + 64] = -n.l[i-1];
    }
  }

  constexpr int8_t& operator[](int l) {
    return next[l + 64];
  }

  constexpr const int8_t& operator[](int l) const {
    return next[l + 64];
  }
};


bool isJoin(const Node& u, const Node& t, int i) {
  NextMap next{t};
  return ((i > 0 && next[u[i-1]] == u[i]) || (i == 0 && t[i] == u[i]));
}

struct ScoredNode {
  Node n;
  int score;
  int id = maxId++;

  static int maxId;
  bool operator<(const ScoredNode& sn) const {
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

  NextMap nextT{t};
  NextMap nextS{s};

  int lastDist = 0;
  int count = 0;
  while (!queue.empty()) {
    ScoredNode us = queue.top();
    Node u = us.n;

    if (count % 1024 == 0 && dist[u] > lastDist) {
      lastDist = dist[u];
      cout << "Finished " << lastDist << ", " << count << ", " << queue.size() << ", " << us.score << endl;
    }
    count++;

    if (u == t) {
      cout << "Finished " << dist[u] << ", " << count <<  ", " << queue.size() << ", " << us.score << endl;
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
      if ((i > 0 && nextT[u[i-1]] == u[i]) || (i == 0 && t[i] == u[i])) {
        continue;
      }

      Node v = u.flip(i);
      if (prev.find(v) != prev.end()) {
        continue;
      }

      // Don't join cut
      if ((i > 0 && nextS[v[i-1]] == v[i]) || (i == 0 && s[i] == v[i])) {
        continue;
      }

      int d = dist[u] + 1;
      if (dist.find(v) == dist.end() || d < dist[v]) {
        prev[v] = i;
        dist[v] = d;
        bool joined = (i > 0 && nextT[v[i-1]] == v[i]) || (i == 0 && t[i] == v[i]);
        bool cut = (i > 0 && nextS[u[i-1]] == u[i]) || (i == 0 && s[i] == u[i]);
        int score = us.score + 1 - joined;
        queue.push({v, score});
      }
    }
  }

  throw runtime_error("no path");
}

vector<int> bfs(Node s, Node t) {
  NextMap smap{s};

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
        bool match = true;
        for (int k = 0; k < u.size(); k++) {
          if (match && u[k] == t[k]) {
            continue;
          }

          match = false;

          // ASSUMES GOAL IS ALPHABETICAL LOWER CASE FOR CORRECTNESS
          // SLOWER IF NOT THE FIRST LETTERS
          if (k > 0 && u[k] == u[k-1] + 1) {
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
        bool match = true;

        for (int k = 0; k < u.size(); k++) {
          if (match && u[k] == s[k]) {
            continue;
          }

          match = false;

          if (k > 0 && smap[u[k-1]] == u[k]) {
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

int main() {
  Node start = parse(input);
  Node fin = start.canonical();
  cout << unparse(start) << ", " << unparse(fin) << endl;
  auto path = astar(start, fin);
  // auto path = bfs(start, fin);

  cout << unparse(start) << endl;
  Node n = start;
  for (int i : path) {
    bool cut = isJoin(n, start, i);
    n = n.flip(i);
    bool join = isJoin(n, fin, i);
    cout << i << ", " << join << ", " <<  cut << endl;
    cout << unparse(n) << endl;
  }
  cout << path.size() << " steps" << endl;
  return 0;
}
