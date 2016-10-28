#include <iostream>
#include <vector>
#include <random>
#include <map>
#include <time.h>
#include <algorithm>

using namespace std;

default_random_engine gen(time(0));
const int TOTAL_CASTLES = 10;
const int TOTAL_SOLDIERS = 100;
array<array<int, TOTAL_SOLDIERS + 1>, TOTAL_CASTLES> wins;

vector<int> random_move() {
    vector<int> x;
    int soldiers = TOTAL_SOLDIERS;
    for (int i = 0; i < TOTAL_CASTLES; i++) {
        discrete_distribution<int> dist(wins[i].begin(), wins[i].begin() + soldiers + 1);
        int s = dist(gen);
        x.push_back(s);
        soldiers -= s;
    }

    if (soldiers > 0) {
        x[TOTAL_CASTLES - 1] += soldiers;
    }

    return x;
}


int play(const vector<int>& a, const vector<int>& b) {
    int score = 0;
    for (int i = 0; i < TOTAL_CASTLES; i++) {
        if (a[i] > b[i]) {
            score += i + 1;
        } else if (b[i] > a[i]) {
            score -= i + 1;
        }
    }

    return score;
}

vector<int> move(int sims, int max_pool) {
    multimap<vector<int>, int> moves;

    for (int i = 0; i < sims; i++) {
        vector<int> move = random_move();
        int total = 0;
        for (auto&& entry : moves) {
            int score = play(move, entry.first);
            if (score > 0) {
                total += 2;
            } else if (score == 0) {
                total += 1;
                entry.second += 1;
            } else {
                entry.second += 2;
            }
        }

        for (int i = 0; i < TOTAL_CASTLES; i++) {
            wins[i][move[i]] += total;
        }
        moves.insert(make_pair(move, total));

        if (moves.size() > max_pool) {
            auto min = *min_element(
                moves.begin(),
                moves.end(),
                [](const auto& p1, const auto& p2) {
                    return p1.second < p2.second;
                }
            );

            moves.erase(min.first);
            for (auto&& entry : moves) {
                int score = play(min.first, entry.first);
                if (score > 0) {
                } else if (score == 0) {
                    entry.second -= 1;
                } else {
                    entry.second -= 2;
                }
            }
        }
    }

//    cout << "{";
//    for (auto& t : moves) {
//        cout << "[";
//        for (auto& s : t.first) {
//            cout << s << ",";
//        }
//        cout << "]";
//        cout << ": " << t.second << ", " << endl;
//    }
//    cout << "}" << endl;
//
//    for (int i = 0; i < TOTAL_CASTLES; i++) {
//        for (int j = 0; j < TOTAL_SOLDIERS; j++) {
//            cout << wins[i][j] << "\t";
//        }
//        cout << endl;
//    }

    return max_element(
        moves.begin(),
        moves.end(),
        [](const auto& p1, const auto& p2) {
            return p1.second < p2.second;
        }
    )->first;
}

int main() {
    int TOTAL_POINTS = TOTAL_CASTLES * (TOTAL_CASTLES + 1) / 2;
    for (int i = 0; i < TOTAL_CASTLES; i++) {
        int x = 1;
        for (int j = 0; j < TOTAL_SOLDIERS + 1; j++) {
            wins[i][j] = max(x, 1);
            if (j * TOTAL_POINTS > TOTAL_CASTLES * TOTAL_SOLDIERS) {
                x--;
            } else {
                x++;
            }
        }
    }

    vector<int> m = move(60000, 50);
    for (auto c : m) {
        cout << c << endl;
    }

    return 0;
}