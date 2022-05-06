#include "bits/stdc++.h"
using namespace std;

typedef long long           lol;
typedef pair<int,int>       pii;
#define pb                  push_back
#define ub                  upper_bound
#define lb                  lower_bound
#define fo(i,l,r,d)         for(auto i=l; d<0?i>r:(d>0?i<r:0); i+=d)
#define all(x)              x.begin(), x.end()
#define ff                  first
#define ss                  second   
#define deb(x)              cout << #x << ' ' << x << '\n'

mt19937 rng (chrono::high_resolution_clock::now().time_since_epoch().count());
template <typename A, typename B> ostream& operator<< (ostream &cout, pair<A, B> const &p) { return cout << "(" << p.first << ", " << p.second << ")"; }
template <typename A, typename B> istream& operator>> (istream& cin, pair<A, B> &p) {cin >> p.first; return cin >> p.second;}
template <typename A> ostream& operator<< (ostream &cout, vector<A> const &v) {cout << "["; for(int i = 0; i < v.size(); i++) {if (i) cout << ", "; cout << v[i];} return cout << "]";}
template <typename A> ostream& operator<< (ostream &cout, set<A> const &v) {cout << "["; for(auto &x: v) cout << x <<", "; return cout << "]";}
template <typename A> istream& operator>> (istream& cin, vector<A> &x){for(int i = 0; i < x.size()-1; i++) cin >> x[i]; return cin >> x[x.size()-1];}
template <typename A, typename B> A amax (A &a, B b){ if (b > a) a = b ; return a; }
template <typename A, typename B> A amin (A &a, B b){ if (b < a) a = b ; return a; }

void split (vector<string>& a, const string& s, const char& d = ' ') {
    a = vector(1, string(0, 'a'));
    for(auto& c: s){
        if(c == d) a.pb("");
        else a.back().pb(c);
    }
}

using Transaction = vector<pii>;

class Dataset {

  public:
    vector<int64_t> TWU;                // TWU values of single items
    size_t N_transactions, N_items;     // Number of transactions and items
    vector<int32_t> item_map;           // mapping of itemID to new hashvalue
    vector<int32_t> item_unmap;         // mapping of new hashvalue to itemID
    vector<Transaction> List;           // List of Transactions
    vector<int32_t> nx;                 // Pointer to next transaction for merging and projection

  public:
    // Reverse Lexicographic order on Transactions
    constexpr static auto Rev_Lex = [](Transaction& x, Transaction& y) {
        const int &nx = x.size(), &ny = y.size(), &n = min(nx, ny);
        for (int i = 1; i <= n; i++)
            if(const int &ix = x[nx-i].ff, &iy = y[ny-i].ff; ix != iy)
                return ix < iy;
        return nx < ny;
    };

    Dataset (const string File) {
        ifstream Fin (File);
        int32_t max_itemID = 0, item_count = 0;

        auto start = chrono::high_resolution_clock::now();
      // Reading the file
        string s;
        while (getline(Fin, s)) {
            vector<string> split_string, item, util; 
            split(split_string, s, ':'), split(item, split_string[0], ' '), split(util, split_string[2], ' ');
            
            assert(item.size() == util.size());
            Transaction T;
            for (int i = 0; i < item.size(); i++)
                if (int it = stoi(item[i]), ut = stoi(util[i]); ut)
                    T.push_back({it, ut}), max_itemID = max(max_itemID, it);
            if (T.size()) List.push_back(T);
        }

      // Sorting items by TWU
        int64_t twu[max_itemID + 1]{0};
        for_each(all(List), [&](Transaction& T){
            // Transaction utility
            int64_t TU = 0;
            for_each(all(T), [&](pii& x){ TU += x.ss; });
            for_each(all(T), [&](pii& x){ twu[x.ff] += TU; });

            // for(auto& [it, ul]: T) TU += ul;
            // for(auto& [it, ul]: T) twu[it] += TU;
        });

        vector<pair<int64_t, int32_t>> order;
        for (int i = 0; i <= max_itemID; i++)
            if(twu[i]) order.push_back({twu[i], i});

        sort(all(order));

        N_items = order.size(), N_transactions = List.size();
        item_map.resize(max_itemID + 1);
        item_unmap.resize(N_items + 1);

      // Hashing items by TWU order
        for (int i = 0; i < N_items; i++)
            item_map[order[i].ss] = i,
            item_unmap[i] = order[i].ss;

      // Transforming Dataset
        for_each(all(List), [&](Transaction& T){
            for_each(all(T), [&](pii& x){
                x.ff = item_map[x.ff];
            });

            // Sort Transactions by TWU order
            sort(all(T));

            // Compressing items
            Transaction X;
            for (int i = 0, q = 1; i < T.size(); X.push_back(T[i]), i += q) {
                while(i+q < T.size() and T[i].ff == T[i+q].ff)
                    T[i].ss += T[i+q].ss, q++;
            }

            T = X;
        });

      // Sorting Dataset
        sort(all(List), Rev_Lex);

        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        cout << "Dataset preparation time: " << duration.count() << "ms\n" << flush;
    }
};

class AlgoTKEH: public Dataset {

  public:
    int64_t min_util = 1;
    map<pair<int32_t, int32_t>, int64_t> EUCS;
    const int K;

    AlgoTKEH (const string& File, const int k): K(k), Dataset(File){

        // {
        //     // debug
        //     cout << "Transactions\n";
        //     for(auto &T: List) cout << T << '\n';
        // }

        auto start = chrono::high_resolution_clock::now();

        RIU();      // Real item utilities strategy
        CUD();      // Co-utility strategy

        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        cout << "min_util raising time: " << duration.count() << "ms\n" << flush;

        // {
        //     // debug
        //     cout << "Transactions\n";
        //     for(auto &T: List) cout << T << '\n';
        // }

        // {
        //     // debug
        //     cout << min_util << '\n';
        // }

        set<int> pri;
        for(int i = 0; i < N_items; i++)
            pri.insert(i);
        vector util(N_transactions, 0ll);

        start = chrono::high_resolution_clock::now();

        search(List, pri, pri, {}, util);

        stop = chrono::high_resolution_clock::now();

        duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        cout << "Search time: " << duration.count() << "ms\n" << flush;

        cout << "Found " << kPatterns.size() << " HUIs\n";

        for(int i = kPatterns.size(); i; i--) {
            auto [ul, itemset] = kPatterns.top();
            ul = -ul;
            for(auto& item: itemset)
                item = item_unmap[item];

            cout << i << " utility: " << ul << '\n';
            cout << itemset << '\n'; 
            kPatterns.pop();
        }
    }

    priority_queue<pair<int64_t, vector<int>>> kPatterns;

    void search(const vector<Transaction> &D,       // projected dataset List-a
                const set<int> &Pri,                // Primary(a)
                const set<int> &Sec,                // Secondary(a)
                const vector<int> &a,               // The current itemset a
                const vector<int64_t>& util){       // utilities of the items in a in U(a)

        for(auto x: Pri) {
            vector b = a; b.pb(x);

            vector<Transaction> b_D;       

            int64_t U_b = 0,
                    LU_b = 0,
                    SU_b = 0;

            map<int, int64_t> lu_bz;

            vector<int64_t> b_util;


            // iterating over transactions
            for(int j = 0; j < D.size(); j++) {
                auto &T = D[j];
                auto tu = util[j];
                int64_t tlu = tu, tsu = tu;

                Transaction b_T; bool found = 0;
                // iterating over items
                for(auto &[it, ul]: T){
                    if(Sec.find(it) != Sec.end() and it != x) {
                        tsu += ul;
                    }
                    if(found)
                        b_T.pb(pair(it, ul));
                    else{
                        if(x == it)
                            found = true,
                            tu += ul, tsu += ul;
                    }
                    tlu += ul;
                }
                if(found) {
                    // computing lu(b U {z}) for all z > x;
                    int64_t tluz = tu;
                    for(auto &[it, ul]: b_T) tluz += ul;
                    for(auto &[it, ul]: b_T) lu_bz[it] += tluz;

                    U_b += tu, LU_b += tlu, SU_b += tsu;
                    if(b_D.empty() or b_D.back() != b_T)
                        b_D.pb(b_T), b_util.pb(tu);
                    else{
                        b_util.back() += tu;
                        for(int i = 0; i < b_T.size(); i++)
                            b_D.back()[i].ss += b_T[i].ss;
                    }
                }
            }

            if(U_b >= min_util) {
                if(kPatterns.size() == K and -kPatterns.top().ff < U_b)
                    kPatterns.push(pair(-U_b, b)), kPatterns.pop(), amax(min_util, -kPatterns.top().ff);
                else if(kPatterns.size() < K)
                    kPatterns.push(pair(-U_b, b));
            }

            if(SU_b < min_util) continue;

            set<int> pri_b, sec_b;
            for(int x: Sec){
                if(lu_bz[x] >= min_util)
                    sec_b.insert(x);
            }

            map<int, int64_t> su_bz;
            for(int j = 0; j < b_D.size(); j++) {
                auto &T = b_D[j];
                int64_t subz = b_util[j];

                for(int i = T.size() - 1; i > -1; i--) {
                    auto &[it, ul] = T[i];
                    if(Sec.count(it))
                        subz += ul, su_bz[it] += subz;
                }
            }
            for( auto y : sec_b ){
                if( su_bz[y] >= min_util ){
                    pri_b.insert(y) ; 
                }
            }

            search( b_D , pri_b , sec_b , b , b_util );
        }
    }

    // Real Item Utilities Strategy
    void RIU () {
        vector U(N_items, int64_t(0));
        for_each(all(List), [&](Transaction& T){
            for_each(all(T), [&](pii& x){
                U[x.ff] += x.ss;
            });
        });
        sort(all(U));
        min_util = max(min_util, U[N_items - min(K, (int)N_items)]);
    }

    // CUD strategy
    map<pair<int, int>, int64_t> CUDM;
    void CUD () {
        for_each(all(List), [&](Transaction &T){
            const int &n = T.size();
            for(int i = 0; i < n; i++)
                for(int j = i+1; j < n; j++)
                    CUDM[pii(T[i].ff, T[j].ff)] += T[i].ss + T[j].ss;
        });

        priority_queue<int64_t> q;
        for(auto &[p, v]: CUDM){
            q.push(-v);
            if(q.size() > K)
                q.pop();
        }

        if(q.size() == K)
            min_util = max(min_util, -q.top());
    }

};


int main(int argc, char** argv){
    ios_base::sync_with_stdio(0), cin.tie(0);

    if(argc < 3) return cout << "Too few arguments\n", 0;

    string file = argv[1];
    const int k = atoi(argv[2]);

    AlgoTKEH(file, k);

    // AlgoTKEH z("example2.txt", 2);
    // AlgoTKEH z("exampleTKEH.txt", 20);
    // AlgoTKEH z("fruithut_utility.txt", 20);

}