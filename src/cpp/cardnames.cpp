#include <algorithm>
#include <array>
#include <cstdint>
#include <deque>
#include <iostream>
#include <queue>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

/**
 * @brief how would you make a solver for this
 * so you get a set of 30 cards out of the 52 in the deck
 * and you have to arrange them in a 5x5 grid, any unused cards are discarded
 * and irrelevant of the ones in the 5x5 grid, the 5 rows and 5 columns are
 * evaluated if a card is used as part of a trick (5-card flush, 5-card
 * straight, pair, trip, quad) it counts as 0 if it's not, it's worth its face
 * value (2-9 = 2-9, 10-K = 10, A = 15) solve for the lowest score standard misc
 * rules: and if a card is not used in both its row and column it counts twice,
 * Aces can be high or low but not wrap (so TJQKA and A2345 are valid but QKA23
 * is not) if the middle 3x3 are 9 distinct values, your final score is halved
 */

/**
 * cards = num * 4 + suit
 * A = 0 and 2-K are 1-12
 */

constexpr std::array<int, 30> CARDS = {1,  3,  4,  5,  7,  9,  11, 12, 15, 17,
                                       18, 20, 21, 23, 25, 26, 27, 29, 32, 35,
                                       36, 38, 39, 41, 42, 45, 46, 47, 48, 49};

int get_suit(int card) { return card & 0x3; }

int get_num(int card) { return card >> 2; }

int get_value(int card) {
  switch (get_num(card)) {
  case 0:
    return 15;
  case 1:
    return 2;
  case 2:
    return 3;
  case 3:
    return 4;
  case 4:
    return 5;
  case 5:
    return 6;
  case 6:
    return 7;
  case 7:
    return 8;
  case 8:
    return 9;
  case 9:
    return 10;
  case 10:
    return 10;
  case 11:
    return 10;
  case 12:
    return 10;
  }
  return 10000;
}

std::string get_str_num(int card) {
  switch (get_num(card)) {
  case 0:
    return "A";
  case 1:
    return "2";
  case 2:
    return "3";
  case 3:
    return "4";
  case 4:
    return "5";
  case 5:
    return "6";
  case 6:
    return "7";
  case 7:
    return "8";
  case 8:
    return "9";
  case 9:
    return "T";
  case 10:
    return "J";
  case 11:
    return "Q";
  case 12:
    return "K";
  }
  return "X";
}

std::string get_str_suit(int card) {
  switch (get_suit(card)) {
  case 0:
    return "C";
  case 1:
    return "D";
  case 2:
    return "H";
  case 3:
    return "S";
  }
  return "Y";
}

std::string get_str(int card) { return get_str_num(card) + get_str_suit(card); }

std::unordered_map<int, std::string> get_card_to_str() {
  std::unordered_map<int, std::string> card_to_str;
  for (int i = 0; i < 52; i++) {
    card_to_str[i] = get_str(i);
  }
  return card_to_str;
}

std::unordered_map<std::string, int> get_str_to_card() {
  std::unordered_map<std::string, int> str_to_card;
  for (int i = 0; i < 52; i++) {
    str_to_card[get_str(i)] = i;
  }
  return str_to_card;
}

const auto CARD_TO_STR = get_card_to_str();
const auto STR_TO_CARD = get_str_to_card();

// precondition: a < b < c < d < e
int get_hand(int a, int b, int c, int d, int e) {
  return (a << 24) + (b << 18) + (c << 12) + (d << 6) + e;
}

using Hand = std::tuple<int, int, int, int, int>;
Hand break_hand(int h) {
  return {h >> 24, (h >> 18) & 63, (h >> 12) & 63, (h >> 6) & 63, h & 63};
}

bool is_flush(int a, int b, int c, int d, int e) {
  int suit = get_suit(a);
  for (int i : {b, c, d, e}) {
    if (suit != get_suit(i)) {
      return false;
    }
  }
  return true;
}

bool is_straight(int a, int b, int c, int d, int e) {
  a = get_num(a);
  b = get_num(b);
  c = get_num(c);
  d = get_num(d);
  e = get_num(e);
  if (b == a + 1 && c == b + 1 && d == c + 1 && e == d + 1) {
    return true;
  }

  if (a != 0) {
    return false;
  }

  return c == b + 1 && d == c + 1 && e == d + 1 && e == 12;
}

int get_value(int a, int b, int c, int d, int e) {
  if (is_flush(a, b, c, d, e)) {
    return 0;
  }

  if (is_straight(a, b, c, d, e)) {
    return 0;
  }

  int score = 0;
  if (get_num(a) != get_num(b)) {
    score += get_value(a);
  }
  if (get_num(b) != get_num(a) && get_num(b) != get_num(c)) {
    score += get_value(b);
  }
  if (get_num(c) != get_num(b) && get_num(c) != get_num(d)) {
    score += get_value(c);
  }
  if (get_num(d) != get_num(c) && get_num(d) != get_num(e)) {
    score += get_value(d);
  }
  if (get_num(e) != get_num(d)) {
    score += get_value(e);
  }
  return score;
}

// precondition: values in x are sorted
std::unordered_map<int, int> get_scores(const int *x, int n) {
  std::unordered_map<int, int> hand_to_score;
  for (int a = 0; a < n; a++) {
    for (int b = a + 1; b < n; b++) {
      for (int c = b + 1; c < n; c++) {
        for (int d = c + 1; d < n; d++) {
          for (int e = d + 1; e < n; e++) {
            int score =
                get_value(*(x + a), *(x + b), *(x + c), *(x + d), *(x + e));
            int hand =
                get_hand(*(x + a), *(x + b), *(x + c), *(x + d), *(x + e));
            hand_to_score.emplace(hand, score);
          }
        }
      }
    }
  }
  return hand_to_score;
}

// map from partial hands to full hands, in order of score (least to greatest)
struct MemMaps {
  std::vector<int> zero;
  std::unordered_map<int, std::vector<int>> one;
  std::unordered_map<int, std::vector<int>> two;
  std::unordered_map<int, std::vector<int>> three;
  std::unordered_map<int, std::vector<int>> four;
  std::unordered_map<int, int> hand_to_score;
};

MemMaps get_mem_maps(const int *x, int n) { return {}; }

/**
 * @brief Solver
 *
 * First do top row, and left col
 * Then do bot row and right col.
 * Then do the middle
 */

constexpr std::array<int, 5> get_null_5() { return {-1, -1, -1, -1, -1}; }

constexpr std::array<std::array<int, 5>, 5> get_null_5by5() {
  std::array<std::array<int, 5>, 5> matrix{};
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      matrix[i][j] = -1;
    }
  }
  return matrix;
}

struct State {
  std::array<int, 5> rows{get_null_5()};
  std::array<int, 5> cols{get_null_5()};
  std::array<std::array<int, 5>, 5> cards{get_null_5by5()};
  int score = 0;
  bool halved = false;
};

int solve(State *state, State *best, int alpha) { return 0; }

/**
 * @brief Here be tests
 *
 * Tests are below!!
 */

void test_gets() {
  for (int i = 0; i < 52; i++) {
    cout << get_suit(i) << ", " << get_num(i) << ", " << get_value(i) << "\n";
  }
}

std::vector<int> get_deck() {
  std::vector<int> v;
  for (int i = 0; i < 52; i++) {
    v.emplace_back(i);
  }
  return v;
}

std::mt19937 get_rand_gen() {
  std::random_device rd;
  std::mt19937 g(rd());
  return g;
}

void test_hands() {
  auto g = get_rand_gen();
  auto v = get_deck();
  for (int i = 0; i < 100; i++) {
    std::shuffle(v.begin(), v.end(), g);
    for (int j = 0; j < 5; j++) {
      cout << v[j] << ", ";
    }
    cout << endl;
    int h = get_hand(v[0], v[1], v[2], v[3], v[4]);
    auto [a, b, c, d, e] = break_hand(h);
    cout << a << ", " << b << ", " << c << ", " << d << ", " << e << endl;
    cout << h << endl;
  }
}

void test_rand_scores() {
  auto hand_to_score = get_scores(CARDS.data(), CARDS.size());
  // 142506
  cout << hand_to_score.size() << endl;

  auto g = get_rand_gen();
  std::vector<int> v{CARDS.begin(), CARDS.end()};
  for (int i = 0; i < 20; i++) {
    std::shuffle(v.begin(), v.end(), g);
  std:
    sort(v.begin(), v.begin() + 5);
    int h = get_hand(v[0], v[1], v[2], v[3], v[4]);
    auto [a, b, c, d, e] = break_hand(h);
    cout << CARD_TO_STR.at(a) << ", ";
    cout << CARD_TO_STR.at(b) << ", ";
    cout << CARD_TO_STR.at(c) << ", ";
    cout << CARD_TO_STR.at(d) << ", ";
    cout << CARD_TO_STR.at(e) << ", ";
    cout << h << ": " << hand_to_score[h] << endl;
  }
}

void test_scores() {
  auto DECK = get_deck();
  // 2598960
  auto hand_to_score = get_scores(DECK.data(), DECK.size());

  cout << hand_to_score.size() << endl;
  std::vector<std::vector<std::string>> v{
      {"3D", "9C", "QD", "QH", "KD"}, {"3S", "6D", "TS", "JD", "QH"},
      {"3C", "3D", "QD", "QH", "KD"}, {"3C", "3D", "QD", "QH", "QS"},
      {"3S", "6S", "TS", "JS", "QS"}, {"3C", "6C", "TC", "JC", "QC"},
      {"4C", "7D", "7S", "QD", "QH"}, {"AC", "7D", "7H", "7S", "QH"},
      {"7H", "8D", "9C", "TH", "QD"}, {"7H", "8D", "9C", "TH", "JD"},
      {"AH", "2D", "3C", "4H", "5D"}, {"AH", "TD", "JC", "QH", "KD"},
      {"AD", "5H", "TS", "JD", "QD"},
  };

  std::vector<int> expected{22, 39, 13, 0, 0, 0, 4, 25, 44, 0, 0, 0, 50};

  for (int i = 0; i < v.size(); i++) {
    std::vector<int> x;
    for (auto &s : v[i]) {
      x.emplace_back(STR_TO_CARD.at(s));
    }
    int h = get_hand(x[0], x[1], x[2], x[3], x[4]);
    cout << h << ": " << hand_to_score[h] << " = " << expected[i] << endl;
  }
}

int main() {
  // test_hands();
  test_rand_scores();
  test_scores();
  // // generate random 30 numbers
  // std::vector<int> v;
  // for (int i = 0; i < 52; i++) {
  //     v.emplace_back(i);
  // }
  // std::random_device rd;
  // std::mt19937 g(rd());

  // std::shuffle(v.begin(), v.end(), g);
  // std::sort(v.begin(), v.begin() + 30);
  // for (int i = 0; i < 30; i++) {
  //     cout << v[i] << ", ";
  // }

  return 0;
}
