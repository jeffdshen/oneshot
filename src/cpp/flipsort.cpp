#include <algorithm>
#include <array>
#include <cstdint>
#include <deque>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

constexpr char input[] = "FlipSortMeZqjkYXWVUDCBghAN";
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

struct node {
  array<int8_t, sizeof(input) - 1> l;

  node flip(int k) const {
    node n;
    for (int i = 0; i < k ; i++) {
      n.l[i] = l[i];
    }

    for (int i = k; i < l.size(); i++) {
      n.l[i] = -l[l.size() - 1 + k - i];
    }

    return n;
  }

  node canonical() const {
    node n;
    for (int i = 0; i < l.size(); i++) {
      n.l[i] = abs(l[i]);
    }

    sort(n.l.begin(), n.l.end());
    return n;
  }

  bool operator==(const node &n) const { return l == n.l; }

  bool operator!=(const node &n) const { return l != n.l; }
};

ostream &operator<<(ostream &os, const node &n) {
  os << "[";
  for (int k : n.l) {
    os << k << ",";
  }
  os << "]";
  return os;
}

namespace std {
template <> struct hash<node> {
  size_t operator()(const node &s) const noexcept {
    size_t h = 0;
    for (int i = 0; i < s.l.size(); i++) {
      h  = ((h << 5) + h) + s.l[i];
    };
    return h;
  }
};
} // namespace std

vector<int> bfs(node s, node t) {
  std::array<int8_t, 256> smap = {0};
  for (int i = 1; i < s.l.size(); i++) {
    smap[s.l[i-1] + 64] = s.l[i];
    smap[-s.l[i] + 64] = -s.l[i-1];
  }

  deque<node> sb;
  deque<node> tb;
  sb.emplace_back(s);
  tb.emplace_back(t);

  if (s == t) {
    return {};
  }

  unordered_map<node, int> sseen;
  sseen.emplace(s, -1);
  unordered_map<node, int> tseen;
  tseen.emplace(t, -1);
  for (int i = 0;; i++) {
    {
      size_t size = sb.size();
      for (int j = 0; j < size; j++) {
        node u = sb.front();
        sb.pop_front();
        bool match = true;
        for (int k = 0; k < u.l.size(); k++) {
          if (match && u.l[k] == t.l[k]) {
            continue;
          }

          match = false;

          // ASSUMES GOAL IS ALPHABETICAL LOWER CASE
          if (k > 0 && u.l[k] == u.l[k-1] + 1) {
            continue;
          }

          node v = u.flip(k);
          if (sseen.find(v) != sseen.end()) {
            continue;
          }

          if (tseen.find(v) != tseen.end()) {
            vector<int> flips;
            flips.resize(i * 2 + 1);
            node p = v;
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
        node u = tb.front();
        tb.pop_front();
        bool match = true;

        for (int k = 0; k < u.l.size(); k++) {
          if (match && u.l[k] == s.l[k]) {
            continue;
          }

          match = false;

          if (k > 0 && smap[u.l[k-1] + 64] == u.l[k]) {
            continue;
          }

          node v = u.flip(k);
          if (tseen.find(v) != tseen.end()) {
            continue;
          }

          if (sseen.find(v) != sseen.end()) {
            vector<int> flips;
            flips.resize(i * 2 + 2);
            node p = v;
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

node parse(const std::string &s) {
  node n;
  for (int i = 0; i < n.l.size(); i++) {
    if (s[i] >= 'a') {
      n.l[i] = s[i] - 'a' + 1;
    } else {
      n.l[i] = -(s[i] - 'A' + 1);
    }
  }
  return n;
}

string unparse(node n) {
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

int main() {
  node start = parse(input);
  // //FlipSortMeZqjkYXWVUDCBghAN
  // //01234567890123456789012345
  // start = start.flip(14);
  // start = start.flip(11);
  // //FlipSortMeZYXWVUDCBghANKJQ
  // //01234567890123456789012345
  // start = start.flip(4);
  // //FlipqjknaHGbcduvwxyzEmTROs
  // //01234567890123456789012345
  // start = start.flip(23);
  // //FlipqjknaHGbcduvwxyzEmTSor
  // //01234567890123456789012345
  // start = start.flip(24);
  // //FlipqjknaHGbcduvwxyzEmTSRO
  // //01234567890123456789012345
  // start = start.flip(8);
  // //FlipqjknorstMeZYXWVUDCBghA
  // //01234567890123456789012345
  // start = start.flip(7);
  // start = start.flip(21);
  // //FlipqjkaHGbcduvwxyzEmnorst
  // //01234567890123456789012345
  // start = start.flip(19);
  // start = start.flip(1);
  // //FEmnorstZYXWVUDCBghAKJQPIL
  // //01234567890123456789012345
  // start = start.flip(24);
  // start = start.flip(22);
  // //FEmnorstZYXWVUDCBghAKJILpq
  // //01234567890123456789012345

  node fin = start.canonical();
  cout << unparse(start) << endl;
  cout << start << ", " << fin << endl;
  auto path = bfs(start, fin);
  node n = start;
  for (int i : path) {
    n = n.flip(i);
    cout << i << endl;
    cout << unparse(n) << endl;
  }
  cout << endl;
  return 0;
}
