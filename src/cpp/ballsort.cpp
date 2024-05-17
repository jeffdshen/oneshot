#include <algorithm>
#include <iostream>
#include <map>
#include <random>
#include <time.h>
#include <vector>

struct Top {
  int open;
  int color;
  int count;
};

using State = std::vector<int>;
using Move = std::pair<int, int>;

struct Solver {
  int width;
  int height;

  std::vector<int> canonical(const std::vector<int>& x) {
    std::vector<std::vector<int>> bins;
    bins.reserve(width);
    for (size_t a = 0; a < width; a++) {
      bins.emplace_back();
      bins.back().reserve(height);
      for (size_t i = 0; i < height; i++) {
        bins.back().emplace_back(x[a * height + i]);
      }
    }
    std::sort(bins.begin(), bins.end());
    std::vector<int> flat;
    flat.reserve(width * height);
    for (auto& bin : bins) {
      for (auto& i : bin) {
        flat.emplace_back(i);
      }
    }
    return flat;
  }

  Top top(const std::vector<int>& x, int a) {
    Top result{};
    for (int i = 1; i <= height; i++) {
      int color = x[(a + 1) * height - i];
      if (color == 0) {
        result.open++;
        continue;
      }

      if (result.color == 0) {
        result.color = color;
        result.count = 1;
        continue;
      }

      if (result.color == color) {
        result.count++;
        continue;
      }

      break;
    }
    return result;
  }

  std::vector<int> make_move(const std::vector<int>& x, int a, int b) {
    auto atop = top(x, a);
    auto btop = top(x, b);
    if (atop.color == 0) {
      return {};
    }

    if (btop.color != 0 && atop.color != btop.color) {
      return {};
    }

    if (atop.count > btop.open) {
      return {};
    }

    std::vector<int> y = x;
    for (int i = atop.count; i >= 1; i--) {
      y[(a + 1) * height - atop.open - i] = 0;
      y[(b + 1) * height - btop.open + atop.count - i] = atop.color;
    }

    return y;
  }

  bool is_solved(const std::vector<int>& x) {
    for (size_t a = 0; a < width; a++) {
      for (size_t i = 0; i < height; i++) {
        if (x[a * height + i] != x[a * height]) {
          return false;
        }
      }
    }
    return true;
  }

  std::vector<Move> get_moves(
      const std::vector<int>& end,
      const std::map<State, std::pair<State, Move>> prev) {
    std::vector<int> v = end;
    std::vector<int> canon_v = canonical(v);
    std::vector<std::pair<int, int>> path;

    std::cout << "reversing" << std::endl;

    while (true) {
      auto it = prev.find(canon_v);
      if (it == prev.end()) {
        break;
      }

      auto [u, move] = it->second;
      path.emplace_back(move);
      v = std::move(u);
      canon_v = canonical(v);
    }
    std::reverse(path.begin(), path.end());
    return path;
  }

  std::vector<std::pair<int, int>> bfs(const std::vector<int>& root) {
    if (is_solved(root)) {
      return {};
    }
    std::vector<State> q;
    q.emplace_back(root);

    using PrevState = std::pair<State, Move>;
    // canonical to prev state
    std::map<State, PrevState> prev;
    int steps = 0;
    while (!q.empty()) {
      std::cout << "next_step " << steps << " " << q.size() << std::endl;
      steps += 1;
      std::vector<State> next_q;
      for (auto& u : q) {
        for (size_t i = 0; i < width; i++) {
          for (size_t j = 0; j < width; j++) {
            if (i == j) {
              continue;
            }
            auto v = make_move(u, i, j);
            if (v.empty()) {
              continue;
            }

            auto canon_v = canonical(v);
            if (prev.find(canon_v) != prev.end()) {
              continue;
            }

            prev.emplace(std::move(canon_v), PrevState{u, {i, j}});
            if (is_solved(v)) {
              return get_moves(v, prev);
            }
            next_q.emplace_back(std::move(v));
          }
        }
      }
      q = std::move(next_q);
    }

    return {};
  }

  void print(std::vector<int>& x) {
    for (size_t i = 0; i < x.size(); i++) {
      if (i % height == 0) {
        std::cout << "| ";
      }
      std::cout << x[i] << ' ';
    }
    std::cout << std::endl;
  }
};

int main() {
  int width;
  int height;
  std::cin >> width;
  std::cin >> height;
  std::vector<int> data;
  data.reserve(width * height);

  for (size_t i = 0; i < width * height; i++) {
    int x;
    std::cin >> x;
    data.emplace_back(x);
  }

  Solver solver{width, height};
  auto path = solver.bfs(data);
  for (auto [a, b] : path) {
    std::cout << (a + 1) << " " << (b + 1) << std::endl;
  }

  return 0;
}