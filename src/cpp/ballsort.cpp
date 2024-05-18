#include <algorithm>
#include <iostream>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <vector>

struct Top {
  int open;
  int color;
  int count;
};

using Bin = std::vector<int>;

struct Slab {
  std::vector<Bin> items;
  std::vector<Top> tops;
  std::vector<int> min_scores;
  std::map<Bin, size_t> index;

  size_t add(Bin item) {
    auto [it, success] = index.emplace(std::move(item), items.size());
    if (success) {
      items.emplace_back(it->first);
      tops.emplace_back(compute_top(it->first));
      min_scores.emplace_back(compute_heuristic(it->first));
    }
    return it->second;
  }

  const Bin& get(size_t i) { return items[i]; }

  const Top& get_top(size_t i) { return tops[i]; }

  int get_min_score(size_t i) { return min_scores[i]; }

  Top compute_top(const Bin& x) {
    Top result{};

    for (size_t i = 1; i <= x.size(); i++) {
      int color = x[x.size() - i];
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

  int compute_heuristic(const Bin& x) {
    int score = 0;
    int prev = 0;
    for (int i = 1; i <= x.size(); i++) {
      auto cur = x[x.size() - i];
      if (cur != prev && prev != 0) {
        score++;
      }
    }
    return score;
  }
};

using State = std::vector<size_t>;
using Move = std::pair<int, int>;

struct BfsState {
  State prev_state;
  Move prev_move;
};

struct AstarState {
  State prev_state;
  Move prev_move;
  int dist;
  bool expanded;
};

struct Solver {
  Slab slab;
  size_t max_solvable = 0;

  State init(std::vector<Bin> x) {
    std::set<int> balls;

    for (auto& bin : x) {
      for (auto& i : bin) {
        balls.emplace(i);
      }
    }
    for (int i : balls) {
      slab.add(Bin(x[0].size(), i));
    }
    max_solvable = balls.size();
    State y;
    y.reserve(x.size());
    for (auto& bin : x) {
      y.emplace_back(slab.add(std::move(bin)));
    }

    return y;
  }

  State canonical(const State& x) {
    State y = x;
    std::sort(y.begin(), y.end());
    return y;
  }

  State make_move(const State& x, int a, int b) {
    auto atop = slab.get_top(x[a]);
    auto btop = slab.get_top(x[b]);
    if (atop.color == 0) {
      return {};
    }

    if (btop.color != 0 && atop.color != btop.color) {
      return {};
    }

    if (atop.count > btop.open) {
      return {};
    }

    Bin ya = slab.get(x[a]);
    Bin yb = slab.get(x[b]);
    for (int i = atop.count; i >= 1; i--) {
      ya[ya.size() - atop.open - i] = 0;
      yb[yb.size() - btop.open + atop.count - i] = atop.color;
    }

    State y = x;
    y[a] = slab.add(std::move(ya));
    y[b] = slab.add(std::move(yb));
    return y;
  }

  std::vector<std::pair<State, Move>> get_moves(const State& u) {
    std::vector<std::pair<State, Move>> results;
    for (size_t i = 0; i < u.size(); i++) {
      for (size_t j = 0; j < u.size(); j++) {
        if (i == j) {
          continue;
        }
        auto v = make_move(u, i, j);
        if (v.empty()) {
          continue;
        }
        results.emplace_back(std::move(v), Move{i, j});
      }
    }
    return results;
  }

  bool is_solved(const State& x) {
    for (auto i : x) {
      if (i >= max_solvable) {
        return false;
      }
    }
    return true;
  }

  template <class SearchState>
  std::vector<Move> get_path(
      const State& end, const std::map<State, SearchState>& prev) {
    State v = end;
    State canon_v = canonical(v);
    std::vector<Move> path;

    // std::cout << "reversing" << std::endl;

    while (true) {
      auto it = prev.find(canon_v);
      if (it == prev.end()) {
        break;
      }
      auto search_state = it->second;
      if (search_state.prev_state.empty()) {
        break;
      }

      path.emplace_back(search_state.prev_move);
      v = std::move(search_state.prev_state);
      canon_v = canonical(v);
    }
    std::reverse(path.begin(), path.end());
    return path;
  }

  std::vector<Move> bfs(const std::vector<size_t>& root) {
    if (is_solved(root)) {
      return {};
    }
    std::vector<State> q;
    q.emplace_back(root);

    // canonical to prev state
    std::map<State, BfsState> prev;
    prev.emplace(root, BfsState{{}, {0, 0}});
    int steps = 0;
    while (!q.empty()) {
      // std::cout << "next_step " << steps << " " << q.size() << std::endl;
      steps += 1;
      std::vector<State> next_q;
      for (auto& u : q) {
        auto vs = get_moves(u);
        for (auto& [v, move] : vs) {
          auto canon_v = canonical(v);
          if (prev.find(canon_v) != prev.end()) {
            continue;
          }

          prev.emplace(std::move(canon_v), BfsState{u, move});
          if (is_solved(v)) {
            return get_path(v, prev);
          }
          next_q.emplace_back(std::move(v));
        }
      }
      q = std::move(next_q);
    }

    return {};
  }

  int get_min_score(const State& x) {
    size_t min_score = 0;
    for (size_t i : x) {
      min_score += slab.get_min_score(i);
    }
    return min_score;
  }

  std::vector<Move> astar(const State& root) {
    if (is_solved(root)) {
      return {};
    }

    // dist + heuristic, state
    using queue_state = std::tuple<int, State>;
    std::priority_queue<
        queue_state,
        std::vector<queue_state>,
        std::greater<queue_state>>
        q;
    std::map<State, AstarState> info;
    info.emplace(canonical(root), AstarState{{}, {0, 0}, 0, false});
    q.emplace(get_min_score(root), root);
    int steps = 0;
    while (!q.empty()) {
      if (steps % 1000 == 0) {
        // std::cout << "next_step " << steps << " " << q.size() << std::endl;
      }
      steps += 1;
      auto next = std::move(q.top());
      q.pop();
      auto& [_, u] = next;
      auto canon_u = canonical(u);
      auto& info_u = info.at(canon_u);
      if (info_u.expanded) {
        continue;
      }
      info_u.expanded = true;
      if (is_solved(u)) {
        return get_path(u, info);
      }

      auto vs = get_moves(u);
      for (auto& [v, move] : vs) {
        int next_dist = info_u.dist + 1;
        auto min_score = get_min_score(v);

        auto canon_v = canonical(v);
        auto [it, success] =
            info.emplace(canon_v, AstarState{{}, move, next_dist + 1, false});
        if (next_dist < it->second.dist) {
          it->second.prev_state = u;
          it->second.prev_move = move;
          it->second.dist = next_dist;
          int f = next_dist + min_score;
          q.emplace(f, v);
        }
      }
    }
    return {};
  }
};

int main() {
  int width;
  int height;
  std::cin >> width;
  std::cin >> height;

  std::vector<Bin> data;
  data.reserve(width);
  for (size_t i = 0; i < width; i++) {
    data.emplace_back();
    data.back().reserve(height);
    for (size_t j = 0; j < height; j++) {
      int x;
      std::cin >> x;
      data.back().emplace_back(x);
    }
  }

  Solver solver;
  auto root = solver.init(data);
  auto path = solver.bfs(root);
  // astar doesn't actually help much since every move just reduces the
  // heuristic by 1.
  // auto path = solver.astar(root);

  for (auto [a, b] : path) {
    std::cout << (a + 1) << " " << (b + 1) << std::endl;
  }

  return 0;
}