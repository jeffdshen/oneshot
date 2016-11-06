/**
 * New version to try to make things more balanced/interesting.
 * Game rules: Each player starts with half the deck, split evenly.
 * Each player puts a card in the pot, which must be strictly higher than the last.
 * A player may also trump the pot with a pair (forcing all others to concede).
 * A player may pass, conceding the pot.
 * The winner of the pot loses all their cards, and takes the rest.
 * The loser must then start a new pot.
 * If a player can no longer make a move, they lose.
 * Generalization: after all passes, the new order is as follows: first to last passer, then next player(s) in order,
 * followed by the winner.
 * Commands: P/p = pass, R/r = resign, x = play card x, xT = play 2 cards X
 *
 **/

#include <map>
#include <unordered_map>
#include <vector>
#include <deque>
#include <set>
#include <iostream>
#include <memory>
#include <math.h>
#include <utility>
#include <random>
#include <algorithm>
#include <fstream>
#include <time.h>

using namespace std;

const bool UNWRITING_ON = true;
const bool LOGGING_ON = true;

class Log : public ostream {
private:
    ofstream fout;
public:
    Log(string s) : fout(s) {}

    template<class T>
    Log& operator<<(const T& x) {
        if (LOGGING_ON) {
            fout << x;
            cout << x;
        }

        return *this;
    }

    typedef std::ostream& (manip)(std::ostream&);
    Log& operator<<(manip& m) {
        if (LOGGING_ON) {
            fout << m;
            cout << m;
        }
        return *this;
    }
};

Log LOG("trumpot.log");

ostream& unwrite(ostream& s, int count) {
    if (LOGGING_ON && UNWRITING_ON) {
        for (int i = 0; i < count; i++) {
            s << "\033[A\33[2K\r" << flush;
        }
    }
    return s;
}


class Move {
private:
    const int card_;
    const int num_;

private:
    Move(int card, int num) : card_(card), num_(num) {}

public:
    int getCard() const {
        return card_;
    }

    int getNum() const {
        return num_;
    }

    bool isPass() const {
        return card_ == 0;
    }

    bool isResign() const {
        return card_ == -1;
    }

    bool isTrump() const {
        return num_ == 2;
    }

    static Move trump(int card) {
        return {card, 2};
    }

    static Move move(int card) {
        return {card, 1};
    }

    static Move pass() {
        return {0, 0};
    }

    static Move resign() {
        return {-1, 0};
    }

    static int encode(const Move& m) {
        return m.card_ * 2 + m.num_ + 2;
    }

    bool operator==(const Move &other) const {
        return card_ == other.card_ && num_ == other.num_;
    }
};

namespace std {
    template<>
    struct hash<Move> {
        size_t operator()(const Move &k) const {
            return Move::encode(k);
        }
    };
}

std::ostream& operator<<(std::ostream &s, const Move &m) {
    if (m.isResign()) {
        return s << "R";
    }

    if (m.isPass()) {
        return s << "P";
    }

    if (m.isTrump()) {
        return s << m.getCard() << "T";
    }

    return s << m.getCard();
}

class Pot {
private:
    vector<vector<int>> pot_;
    int size_;
    int top_;

public:
    Pot(int players) : pot_(players), size_(0) {}

    const vector<vector<int>>& getPot() const {
        return pot_;
    };

    bool empty() const {
        return size_ == 0;
    }

    int top() const {
        return empty() ? 0 : top_;
    };

    void add(int player, int i) {
        pot_[player].push_back(i);
        ++size_;
        top_ = i;
    };

    int getSize(int player) const {
        return pot_[player].size();
    }

    const int getSize() const {
        return size_;
    }

    void remove(int player) {
        pot_[player].clear();
    }

    void clear() {
        for (auto& p : pot_) {
            p.clear();
        }
        top_ = 0;
        size_ = 0;
    }
};


std::ostream& operator<<(std::ostream &s, const Pot &p) {
    for (int i = 0; i < p.getPot().size(); i++) {
        s << "Pot " << i << ":";
        for (auto card : p.getPot()[i]) {
            s << card << " ";
        }

        s << endl;
    }

    return s;
}


class Hand {
private:
    map<int, int> hand_;

public:
    Hand(map<int, int> hand) : hand_(hand) {}

    const map<int, int>& getHand() const {
        return hand_;
    };

    bool empty() const {
        return hand_.empty();
    }

    void add(Pot& pot) {
        for (auto& pp : pot.getPot()) {
            for (auto card : pp) {
                ++hand_[card];
            }
        }
    }

    void remove(Move& m) {
        if (m.isPass()) {
            return;
        }

        int card = m.getCard();
        hand_[card] -= m.getNum();
        if (hand_[card] == 0) {
            hand_.erase(card);
        }
    }

    int top() const {
        return empty() ? 0 : hand_.crbegin()->first;
    }
};

std::ostream& operator<<(std::ostream &s, const Hand &h) {
    int count = 0;
    for (auto ent : h.getHand()) {
        count += ent.second;
        for (int i = 0; i < ent.second; i++) {
            s << ent.first << " ";
        }
    }

    s << "(" << count << ")";

    return s;
}

class Board {
private:
    vector<Hand> hands_;
    set<int> losers_;
    deque<int> order_;
    deque<int> passed_;
    Pot pot_;
    int players_;

private:
    void nextRound() {
        int winner = order_.back();
        pot_.remove(winner);
        hands_[winner].add(pot_);
        pot_.clear();

        while (!passed_.empty()) {
            order_.push_front(passed_.back());
            passed_.pop_back();
        }
    }

public:
    Board(vector<Hand> h) : hands_(h), losers_(), pot_(h.size()), order_(), passed_(), players_(h.size()){
        for (int i = 0; i < h.size(); i++) {
            order_.push_back(i);
        }
    }

    Board(const Board& b) :
            hands_(b.getHands()),
            losers_(b.getLosers()),
            order_(b.getOrder()),
            passed_(b.getPassed()),
            pot_(b.getPot()),
            players_(b.getPlayers()){}

    const vector<Hand>& getHands() const {
        return hands_;
    }

    const set<int>& getLosers() const {
        return losers_;
    }

    const deque<int>& getOrder() const {
        return order_;
    }

    const deque<int>& getPassed() const {
        return passed_;
    }

    const Pot getPot() const {
        return pot_;
    }

    int getPlayers() const {
        return players_;
    }

    int getCurrentPlayer() const {
        return order_.front();
    }

    bool hasLost() const {
        return pot_.empty() && hands_[getCurrentPlayer()].empty();
    }

    bool isGameOver() const {
        return order_.size() == 1 && passed_.empty();
    }

    bool canMove(Move m) const {
        if (m.isResign()) {
            return true;
        }

        if (hasLost()) {
            return false;
        }

        if (m.isPass()) {
            return !pot_.empty();
        }

        if (m.getCard() <= pot_.top()) {
            return false;
        }

        int player = getCurrentPlayer();
        auto& h = hands_[player].getHand();
        if (h.find(m.getCard()) == h.end() || m.getNum() > h.at(m.getCard())) {
            return false;
        }

        return true;
    }

    vector<Move> getMoves() const {
        vector<Move> result;
        result.push_back(Move::resign());
        if (hasLost()) {
            return result;
        }

        if (!pot_.empty()) {
            result.push_back(Move::pass());
        }
        int player = getCurrentPlayer();
        int top = pot_.top();
        for (auto& h : hands_[player].getHand()) {
            if (h.first > top) {
                result.push_back(Move::move(h.first));
                if (h.second >= 2) {
                    result.push_back(Move::trump(h.first));
                }
            }
        }

        return result;
    }

    /**
     * Makes the move on the board by the player
     */
    void move(Move m) {
        int player = getCurrentPlayer();
        if (m.isResign()) {
            losers_.insert(player);
            order_.pop_front();
            if (order_.size() == 1) {
                nextRound();
            }
            return;
        }

        if (m.isPass()) {
            order_.pop_front();
            passed_.push_back(player);
            if (order_.size() == 1) {
                nextRound();
            }

            return;
        }

        hands_[player].remove(m);
        pot_.add(player, m.getCard());

        if (m.isTrump()) {
            order_.pop_front();
            order_.push_back(player);
            nextRound();
        } else {
            order_.pop_front();
            order_.push_back(player);
        }
    };
};


std::ostream& operator<<(std::ostream& s, const Board& b) {
    for (auto i : b.getPassed()) {
        s << "Passed " << i << ": " << b.getHands()[i] << endl;
    }

    for (auto i : b.getOrder()) {
        s << "In " << i << ": " << b.getHands()[i] << endl;
    }

    s << b.getPot();

    return s;
}

std::ostream& unwrite(std::ostream& s, const Board& b) {
    int count = 0;
    for (auto i : b.getPassed()) {
        count++;
    }

    for (auto i : b.getOrder()) {
        count++;
    }

    count += b.getHands().size();
    return unwrite(s, count);
}

class Player {
public:
    virtual ostream& comment(ostream& s) const=0;
    virtual Move move(const Board& board)=0;
};

std::ostream& operator<<(std::ostream& s, const Player& p) {
    return p.comment(s);
}

class CommandLinePlayer : public Player{
public:
    Move move(const Board& board) {
        LOG << "Please enter your move: " << endl;

        for (int count = 1; ; count++) {
            string s;
            getline(cin, s);
            unwrite(cout, 2);
            if (s == "P" || s == "p") {
                return Move::pass();
            }

            if (s == "R" || s == "r") {
                return Move::resign();
            }

            if (s.back() == 'T' || s.back() == 't') {
                s.pop_back();
                try {
                    return Move::trump(stoi(s));
                } catch (exception e) {

                }
            } else {
                try {
                    return Move::move(stoi(s));
                } catch (exception e) {

                }
            }

            LOG << s << " is not a valid move (P/R/x/xT):" << endl;
        }
    }

    ostream& comment(ostream& s) const {
        return s;
    }
};

class Game {
private:
    vector<shared_ptr<Player>> players_;
    Board board_;
public:
    Game(
            vector<shared_ptr<Player>> players, Hand h
    ) : players_(players), board_({players.size(), h}) {}

    Game(
            vector<shared_ptr<Player>> players, vector<Hand> h
    ) : players_(players), board_(h) {}

    int run() {
        while (true) {
            int player = board_.getCurrentPlayer();
            LOG << board_ << endl;
            if (board_.isGameOver()) {
                LOG << "Player " << player << " wins!" << endl;
                return player;
            }

            if (board_.hasLost()) {
                unwrite(cout, 1);
                unwrite(cout, board_);
                LOG << "Player " << player << " lost" << endl;

                board_.move(Move::resign());
                continue;
            }

            LOG << "Player " << player << " is moving..." << endl;

            Move m = players_[player]->move(board_);

            unwrite(cout, 2);
            unwrite(cout, board_);

            if (board_.canMove(m)) {
                board_.move(m);

                LOG << "Player " << player << " played " << m << endl;
                LOG << *players_[player];
            } else {
                LOG << "Player " << player << " made an invalid move: " << m << endl;
            }
        }
    }
};

class RandomPlayout {
private:
    static default_random_engine gen_;
    static discrete_distribution<int> dist_;
    static discrete_distribution<int> beat_;
    static discrete_distribution<int> high_;

public:
    static Move move(const Board& board);

    static int simulate(Board& board) {
        while (!board.isGameOver()) {
            Move m = move(board);
            board.move(m);
        }

        return board.getCurrentPlayer();
    }
};


default_random_engine RandomPlayout::gen_(1337);
discrete_distribution<int> RandomPlayout::dist_{1, 1, 7, 2};
discrete_distribution<int> RandomPlayout::beat_{0, 3, 0, 2};
discrete_distribution<int> RandomPlayout::high_{1, 1, 2, 8};


Move RandomPlayout::move(const Board& board) {
    if (board.hasLost()) {
        return Move::resign();
    }

    auto& hands = board.getHands();
    int player = board.getCurrentPlayer();
    auto& hand = hands[player].getHand();
    auto& h = hands[player];
    auto& pot = board.getPot();

    if (hand.empty()) {
        return Move::pass();
    }

    // if empty, play smallest card
    if (pot.empty()) {
        return Move::move(hand.begin()->first);
    }

    // if top is too big, pass
    if (pot.top() >= h.top()) {
        return Move::pass();
    }

    // play randomly
    int num = 0;
    int players = board.getPlayers();

    // if would get more cards back than lost, trump
    int lost = 2 * pot.getSize(player) - pot.getSize();
    int oppLost = (pot.getSize() - pot.getSize(player)) / (players - 1);
    if (lost + 3 < oppLost) {
        num = 1;
    } else if (lost + 2 < oppLost) {
        num = beat_(gen_);
    } else if (lost + 1 < oppLost) {
        num = high_(gen_);
    } else {
        num = dist_(gen_);
    }

    switch (num) {
        case 0: {
            // pass
            return Move::pass();
        }
        case 1: {
            // play trump
            for (auto it = hand.upper_bound(pot.top()); it != hand.end(); ++it) {
                if (it->second >= 2) {
                    return Move::trump(it->first);
                }
            }

            return Move::pass();
        }
        case 2: {
            // play higher card
            uniform_int_distribution<int> uniform(pot.top() + 1, h.top());
            int pick = uniform(gen_);
            return Move::move(hand.lower_bound(pick)->first);
        }
        case 3: {
            // play highest
            bool winning = true;
            for (int i = 0; i < players; i++) {
                if (hands[i].top() > h.top()) {
                    winning = false;
                    break;
                }
            }

            if (winning) {
                return Move::move(h.top());
            } else {
                return Move::pass();
            }
        }
        default: {
            return Move::pass();
        }
    }
}

class MonteCarloNode {
private:
    static const int kN = 1; // visit 5 times before expanding
    static default_random_engine gen_;
    static uniform_int_distribution<int> dist_;

    int w_; // wins for last player
    int n_; // number of times visited
    int p_; // last player (init by parent)
    unordered_map<Move, shared_ptr<MonteCarloNode>> next_;

public:
    MonteCarloNode() : w_(0), n_(0), p_(0), next_(){}

    double uct(int t);

    Move move() {
        auto element = max_element(next_.begin(), next_.end(), [](const auto& a, const auto& b){
            return a.second->n_ < b.second->n_;
        });

        return element->first;
    }

    double winPercentage(Move& m) {
        auto& x = *next_[m];
        return x.w_ * 1.0 / x.n_;
    }

    void initNext(const Board& board) {
        if (board.isGameOver()) {
            return;
        }

        vector<Move> moves = board.getMoves();
        for (Move& move : moves) {
            next_[move] = make_shared<MonteCarloNode>();
            next_[move]->p_ = board.getCurrentPlayer();
        }

    }

    int select(Board& board);
};

default_random_engine MonteCarloNode::gen_(1776);
uniform_int_distribution<int> MonteCarloNode::dist_(0, 10000);

double MonteCarloNode::uct(int t) {
    if (n_ == 0) {
        return 10000 + dist_(gen_);
    }

    return ((double) w_) / n_ + 0.45 * sqrt(log(t) / n_);
}

int MonteCarloNode::select(Board& board) {
    if (n_ < kN) {
        int winner = RandomPlayout::simulate(board);
        if (p_ == winner) {
            w_++;
        }

        n_++;
        return winner;
    }

    if (n_ == kN) {
        initNext(board);
    }

    if (next_.empty()) {
        int winner = board.getCurrentPlayer();
        if (p_ == winner) {
            w_++;
        }

        n_++;
        return winner;
    }

    auto element = max_element(next_.begin(), next_.end(), [this](auto& a, auto& b){
        return a.second->uct(n_) < b.second->uct(n_);
    });


    board.move(element->first);

    int winner = element->second->select(board);

    if (winner == p_) {
        w_++;
    }
    n_++;
    return winner;
}


class MonteCarloPlayer : public Player {
private:
    const int sims_;
    const int timeLimit_;
    double winP_;

public:
    MonteCarloPlayer(int sims, int timeLimit) : sims_(sims), timeLimit_(timeLimit) {}

    Move move(const Board& board) {
        auto begin = time(0);
        MonteCarloNode root;
         LOG << 0 << "/" << timeLimit_ << endl;

        for (int i = 0; i < sims_; i++) {
            if ((i + 1) % 10000 == 0) {
                auto end = time(0);
                if (end - begin >= timeLimit_) {
                    break;
                }
                unwrite(cout, 1);
                LOG << (end - begin) << "/" << timeLimit_ << " (" << i << "/" << sims_ << ")" << endl;
            }
            Board boardCopy(board);
            root.select(boardCopy);
        }

        unwrite(cout, 1);
        Move m = root.move();
        winP_ = root.winPercentage(m);

        return m;
    }

    ostream& comment(ostream& s) const {
        return s << "Winning probability: " << winP_ << endl;
    }
};

class SimulationPlayer : public Player {
public:
    Move move(const Board& board) {
        return RandomPlayout::move(board);
    }

    ostream& comment(ostream& s) const {
        return s;
    }
};


int main() {
    auto p = make_shared<CommandLinePlayer>();
    auto q = make_shared<MonteCarloPlayer>(1000000, 5);
    auto r = make_shared<SimulationPlayer>();
    vector<shared_ptr<Player>> players{shared_ptr<Player>(q), shared_ptr<Player>(r)};
//    Hand h {{{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 2}, {8, 2}, {9, 2}, {10, 2}, {11, 2}, {12, 2}, {13, 2}}};
    Hand h {{{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 2}, {8, 2}, {9, 2}, {10, 2}}};

    Game g{players, h};

    g.run();
}