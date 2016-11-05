/**
 * Game rules: Each player starts with half the deck, split evenly.
 * Each player puts a card in the pot, which must be strictly higher than the last.
 * A player may also trump the pot with a pair that is at least as high as the last.
 * A player may pass, conceding the pot.
 * The winner of the pot loses all the highest cards in the pot, and takes the rest.
 * The loser must then start a new pot.
 * If a player can no longer make a move, they lose.
 * Generalization: after all passes, the new order is as follows: first to last passer, then next player(s) in order,
 * followed by the winner.
 * Commands: P/p = pass, R/r = resign, x = play card x, xT = play 2 cards X
 *
 * Conclusion: seems that strategy is not too deep - it's too easy to trade 1-1 so,
 * if your opponent ever gains a card advantage, they can just maintain it til the end.
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
    map<int, int> pot_;

public:
    const map<int, int>& getPot() const {
        return pot_;
    };

    bool empty() const {
        return pot_.empty();
    }

    int top() const {
        return empty() ? 0 : pot_.crbegin()->first;
    };

    void add(int i) {
        ++pot_[i];
    };

    void remove() {
        pot_.erase(top());
    }

    void clear() {
        pot_.clear();
    }

    void add(Move& m) {
        if (m.isPass()) {
            return;
        }

        pot_[m.getCard()] += m.getNum();
    }
};


std::ostream& operator<<(std::ostream &s, const Pot &p) {
    for (auto ent : p.getPot()) {
        for (int i = 0; i < ent.second; i++) {
            s << ent.first << " ";
        }
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
        for (auto& entry : pot.getPot()) {
            hand_[entry.first] += entry.second;
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

private:
    void nextRound() {
        int winner = order_.back();
        pot_.remove();
        hands_[winner].add(pot_);
        pot_.clear();

        while (!passed_.empty()) {
            order_.push_front(passed_.back());
            passed_.pop_back();
        }
    }

public:
    Board(vector<Hand> h) : hands_(h), losers_(), pot_(), order_(), passed_() {
        for (int i = 0; i < h.size(); i++) {
            order_.push_back(i);
        }
    }

    Board(const Board& b) :
        hands_(b.getHands()),
        losers_(b.getLosers()),
        order_(b.getOrder()),
        passed_(b.getPassed()),
        pot_(b.getPot()) {}

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

        if (m.getCard() < pot_.top() || (!m.isTrump() &&  m.getCard() == pot_.top())) {
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
            }

            if (h.first >= top && h.second >= 2) {
                result.push_back(Move::trump(h.first));
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
        pot_.add(m);

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

    s << "Pot: " << b.getPot() << endl;

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

    count++;
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
discrete_distribution<int> RandomPlayout::dist_{3, 1, 4};


Move RandomPlayout::move(const Board& board) {
    if (board.hasLost()) {
        return Move::resign();
    }

    auto& hand = board.getHands()[board.getCurrentPlayer()].getHand();
    auto& pot = board.getPot();
    // if empty, play smallest card
    if (pot.empty()) {
        return Move::move(hand.begin()->first);
    }

    // otherwise, either pass, trump with smallest possible, or play random higher card
    int num = dist_(gen_);
    switch (num) {
        case 0: {
            return Move::pass();
        }
        case 1: {
            for (auto it = hand.lower_bound(pot.top()); it != hand.end(); ++it) {
                if (it->second >= 2) {
                    return Move::trump(it->first);
                }
            }

            return Move::pass();
        }
        case 2: {
            uniform_int_distribution<int> uniform(0, hand.size() - 1);
            int pick = uniform(gen_);
            for (auto& entry : hand) {
                if (pick-- == 0) {
                    if (pot.top() <= entry.first) {
                        return Move::move(entry.first);
                    } else {
                        return Move::pass();
                    }
                }
            }

            return Move::pass();
        }
        default: {
            return Move::pass();
        }
    }
}

class MonteCarloNode {
private:
    static const int kN = 5; // visit 5 times before expanding
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
    double winP_;

public:
    MonteCarloPlayer(int sims) : sims_(sims) {}

    Move move(const Board& board) {
        MonteCarloNode root;
        LOG << 0 << "/" << sims_/ 10000 << endl;
        for (int i = 0; i < sims_; i++) {
            if ((i + 1) % 10000 == 0) {
                unwrite(cout, 1);
                LOG << (i + 1) / 10000 << "/" << sims_/ 10000 << endl;
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

int main() {
    auto p = make_shared<CommandLinePlayer>();
    auto q = make_shared<MonteCarloPlayer>(100000);
    vector<shared_ptr<Player>> players{shared_ptr<Player>(q), shared_ptr<Player>(q)};
    Hand h {{{1, 2}, {2, 2}, {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 2}, {8, 2}, {9, 2}, {10, 2}, {11, 2}, {12, 2}, {13, 2}}};
    Game g{players, h};

    g.run();
}