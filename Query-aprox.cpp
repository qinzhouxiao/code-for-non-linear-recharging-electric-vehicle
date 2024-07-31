#include <iostream>
#include <vector>
#include <queue>
#include <time.h>
#include <map>
#include <unordered_map>
#include <cmath>
#include <fstream>
#include <assert.h>
#include <limits.h>
using namespace std;
const int MAXN = 2800000;
const int MAXM = 51 + 2;
const int INF = 0x3f3f3f3f;

struct EDGE {
    int u, v, t, d;
};
EDGE x;

struct LABEL {
    LABEL() {}
    LABEL(int a1, int a2, int a3, int a4) : sum_t(a1), sum_d(a2), pos(a3), index(a4) { }
    int sum_t, sum_d, pos, index;
    friend bool operator >(LABEL a, LABEL b) {
        if (a.sum_d != b.sum_d) return a.sum_d > b.sum_d;
        if (a.sum_t != b.sum_t) return a.sum_t > b.sum_t;
        if (a.pos != b.pos) return a.pos > b.pos;
        return a.index > b.index;
    }
    friend bool operator <(LABEL a, LABEL b) {
        if (a.sum_d != b.sum_d) return a.sum_d < b.sum_d;
        if (a.sum_t != b.sum_t) return a.sum_t < b.sum_t;
        if (a.pos != b.pos) return a.pos < b.pos;
        return a.index < b.index;
    }
    friend bool operator ==(LABEL a, LABEL b) {
        return a.sum_d == b.sum_d && a.sum_t == a.sum_t && a.pos == b.pos && a.index == b.index;
    }
};

int Capacity, initial_d;

void Init() {
    Capacity = 2000000;
    // initial_d = 1000000;
}

int n, m;
int station_num;
int SS, TT;
bool SS_is_station = true, TT_is_station = true;
// vector<EDGE> graph[MAXN];
// vector<EDGE> sketch[MAXN];
vector<vector<EDGE>> graph(MAXN, vector<EDGE>());
vector<vector<EDGE>> sketch(MAXN, vector<EDGE>());
vector<int> stations;
unordered_map<int, int> trans;
unordered_map<int, int> is_station;

// NY BAY COL FLA
// NW NE CAL LKS
ifstream road_data("dataset/graph/USA-road-t-d.NY.gr");
ifstream sketch_data("dataset/NY-road-station.txt");;

ofstream log_out("log/aprox/test-log-NY-Q15.txt");
ofstream exp_out("log/aprox/test-aprox-NY-Q15.txt");


void ReadGraph() {
    std::ios::sync_with_stdio(false);
    
    road_data >> n >> m;
    for (int i = 1;i <= m; i++) {
        EDGE e;
        road_data >> e.u >> e.v >> e.t >> e.d;
        graph[e.u].push_back(e);
    }
    
    
    sketch_data >> station_num;
    for (int i = 1;i <= station_num; i++) {
        int x;
        sketch_data >> x;
        stations.push_back(x);
        trans[x] = i;
        is_station[x] = 1;
    }
    
    for (int i = 1; i <= station_num * (station_num - 1); i++) {
        EDGE e;
        sketch_data >> e.u >> e.v >> e.t >> e.d;
        if (e.t == -1 || e.d == -1) continue;
        sketch[e.u].push_back(e);
    }
    
    // 占位
    stations.push_back(0);
    stations.push_back(0);
    // task_data >> SS >> TT;
    // if (trans[SS] == 0) { stations.push_back(SS); trans[SS] = station_num + 1; SS_is_station = false; }
    // if (trans[TT] == 0) { stations.push_back(TT); trans[TT] = station_num + 2; TT_is_station = false; }
}

int F(int t, int d) {
    assert(d <= Capacity);
    // 0.5h=50%=50km=50000m=1800s
    // 0.5h=20%=20km=20000m=1800s
    //   3h=20%=20km=20000m=10800s
    //   4h=10%=10km=10000m=14400s
    int ret = t;
    ret += min(max(0, d - (int)(0.0 * Capacity)), (int)(0.5 * Capacity)) * 0.36;
    ret += min(max(0, d - (int)(0.5 * Capacity)), (int)(0.2 * Capacity)) * 0.9;
    ret += min(max(0, d - (int)(0.7 * Capacity)), (int)(0.2 * Capacity)) * 5.4;
    ret += min(max(0, d - (int)(0.9 * Capacity)), (int)(0.1 * Capacity)) * 14.4;
    // ret += max(0, d - Capacity) * 1000;
    return ret;
}

// tmin表示以t作为关键字搜索最短路径 0/1分别表示这条路径的时间、路程
// dmin表示以d作为关键字搜索最短路径 0/1分别表示这条路径的时间、路程
// int tmin[MAXM][MAXN][2];
// int dmin[MAXM][MAXN][2];
//1 vector<vector<vector<int>>> tmin(MAXM, vector<vector<int>>(MAXN, vector<int>(2, 0)));
//1 vector<vector<vector<int>>> dmin(MAXM, vector<vector<int>>(MAXN, vector<int>(2, 0)));
unordered_map<int, int> tmin[MAXM][2], dmin[MAXM][2];

int weight_t(EDGE e) { return e.t; }
int weight_d(EDGE e) { return e.d; }

// int dis[MAXN];
vector<int> dis(MAXN, INF);
vector<int> dis_t(MAXN, INF);
vector<int> dis_d(MAXN, INF);
void Dijkstra(int S, unordered_map<int, int> ret[], int (*weight)(EDGE)) {
    
    // clock_t start = clock();
    
    // memset(dis, 0x3f, sizeof(dis));
    fill(dis.begin(), dis.end(), INF);
    fill(dis_t.begin(), dis_t.end(), INF);
    fill(dis_d.begin(), dis_d.end(), INF);
    ret[0].clear();
    ret[1].clear();
    dis[S] = dis_t[S] = dis_d[S] = 0;
    
    priority_queue<pair<int,int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.push(make_pair(dis[S], S));
    
    while (!pq.empty()) {
        pair<int, int> pr = pq.top(); pq.pop();
        if (pr.first != dis[pr.second]) continue;
        
        for (auto e : graph[pr.second]) {
            // if (ret[e.u][1] + e.d >= Capacity) continue;
            int w = weight(e);
            if (dis[e.u] + w < dis[e.v]) {
                dis[e.v] = dis[e.u] + w;
                pq.push(make_pair(dis[e.v], e.v));
                
                dis_t[e.v] = dis_t[e.u] + e.t;
                dis_d[e.v] = dis_d[e.u] + e.d;
            }
        }
    }
    
    for (int i = 1;i <= n; i++) {
        if (dis_d[i] <= Capacity) {
            ret[0][i] = dis_t[i];
            ret[1][i] = dis_d[i];
        }
    }
    
    // clock_t end = clock();

    // log_out << " Dijkstra Time : " << (double)(end - start) / CLOCKS_PER_SEC << " s  " << endl;
}

void PreProcess(int L, int R) {
    
    // log_out << "------   PreProcess Part   ----" << endl;
    // for (int i = 0;i < stations.size(); i++) {
    for (int i = L;i <= R; i++) {
        // log_out << "Processing : " << i << " " << stations[i] << endl;
        Dijkstra(stations[i], tmin[i+1], weight_t);
        Dijkstra(stations[i], dmin[i+1], weight_d);
    }
}

// vector<LABEL> skyline[2][MAXN];
vector<vector<vector<LABEL>>> skyline(2, vector<vector<LABEL>>(MAXN, vector<LABEL>()));
queue<pair<int, int>> dirt;

LABEL ans_label1, ans_label2;
int ans_t, ans_d;
int SNUM;
int flag = 0;
int T_limit, D_limit;

int PQNUM = 0;
int AVLTime = 0;
int TotTime;
int TotPair;
// int Tmin[2][MAXN];
// int Lmin[2][MAXN];
vector<vector<int>> Tmin(2, vector<int>(MAXN, 0));
vector<vector<int>> Lmin(2, vector<int>(MAXN, 0));
int ID[2];

int ans;
int total_ans = INT_MAX;

int FindIndex(LABEL &label) {
    for (int i = 0;i < skyline[label.index][label.pos].size(); i++) {
        if (label == skyline[label.index][label.pos][i]) {
            return i;
        }
    }
    return -1;
}

bool PushLabel(LABEL &label) {
    for (auto &that_label : skyline[label.index][label.pos]) {
        if (label.sum_t >= that_label.sum_t && label.sum_d >= that_label.sum_d) {
            return false;
        }
    }
    
    if (skyline[label.index][label.pos].size() < 1) {
        dirt.push(make_pair(label.index, label.pos));
        skyline[label.index][label.pos].push_back(label);
        return true;
    }
    if (label.sum_t < skyline[label.index][label.pos][0].sum_t) {
        skyline[label.index][label.pos][0] = label;
        return true;
    }
    
    if (skyline[label.index][label.pos].size() < 2) {
        skyline[label.index][label.pos].push_back(label);
        return true;
    }
    if (label.sum_d < skyline[label.index][label.pos][1].sum_d) {
        skyline[label.index][label.pos][1] = label;
        return true;
    }
    
    if (skyline[label.index][label.pos].size() < 15) {
        skyline[label.index][label.pos].push_back(label);
        return true;
    }
    if ((label.sum_t + label.sum_d) % 2) {
        skyline[label.index][label.pos][1] = label;
        return true;
    }
    
    return false;
}

bool LabelCoupling2(LABEL &label, int &ans) {
    if (skyline[!label.index][label.pos].size()) {
        int tmin = skyline[!label.index][label.pos].back().sum_d;
        if (label.sum_t + tmin > T_limit) return true;
        
        for (auto that_label : skyline[!label.index][label.pos]) {
            SNUM++;
            
            if (label.sum_t + that_label.sum_t > T_limit) continue;
            if (label.sum_d + that_label.sum_d > D_limit) break;
//            if (a../ns <= F(label.sum_t + tmin, label.sum_d + that_label.sum_d)) break;
            
            auto subans = F(label.sum_t + that_label.sum_t, label.sum_d + that_label.sum_d);
            if (subans < ans) {
                ans = subans;
                ans_label1 = label;
                ans_label2 = that_label;
                ans_t = label.sum_t + that_label.sum_t;
                ans_d = label.sum_d + that_label.sum_d;
                // log_out << FindIndex(label) << " " << FindIndex(that_label) << endl;
            }
        }
        return true;
    }
    return false;
}

bool LabelPruning2(LABEL &label, int &ans) {
    int dmin_t, dmin_d, tmin_t, tmin_d;
    
    // 到T的最省电的路线都已经超过Capacity
    if (dmin[ID[!label.index]][1].count(label.pos) == 0) {
        return true;
    }   else {
        dmin_t = dmin[ID[!label.index]][0][label.pos];
        dmin_d = dmin[ID[!label.index]][1][label.pos];
    }
    
    // 如果到T的最短路线超过了Capacity，则无法使用这条道路进行pruning
    if (tmin[ID[!label.index]][1].count(label.pos) == 0) {
        tmin_t = tmin_d = 0;
    }   else {
        tmin_t = tmin[ID[!label.index]][0][label.pos];
        tmin_d = tmin[ID[!label.index]][1][label.pos];
    }
    
    
    // 到终点的距离大于Capacity
    if (label.sum_t + tmin_t > T_limit) return true;
    if (label.sum_d + dmin_d > D_limit) return true;

    // 预估的F值已经超过了ans
    if (F(label.sum_t + tmin_t,
          label.sum_d + dmin_d) >= ans) return true;
    
    // 当前的t/d已经大于贪心上界
    if (dmin[ID[label.index]][0].count(label.pos) == 1) {
        if (label.sum_t > dmin[ID[label.index]][0][label.pos]) return true;
    }
    if (tmin[ID[label.index]][1].count(label.pos) == 1) {
        if (label.sum_d > tmin[ID[label.index]][1][label.pos]) return true;
    }
    
    if (label.sum_t > Tmin[label.index][label.pos] && Tmin[label.index][label.pos] != -1) return true;
    
    int newL = F(label.sum_t, label.sum_d + dmin_d);
    if (newL > Lmin[label.index][label.pos] && Lmin[label.index][label.pos] != -1) return true;
    
    return false;
}

void Solve2(int S, int T) {
    int init_sol1, init_sol2;
    if (tmin[trans[S]][1].count(T) == 0) init_sol1 = 0; else init_sol1 = F(tmin[trans[S]][0][T], tmin[trans[S]][1][T]);
    if (dmin[trans[S]][1].count(T) == 0) init_sol2 = 0; else init_sol2 = F(dmin[trans[S]][0][T], dmin[trans[S]][1][T]);
    
    ans = ans_t = ans_d = -1;
    
    if (init_sol2 == 0) {
        TotPair--;
        return;
        
    }   else if (init_sol1 == 0) {
        T_limit = dmin[trans[S]][0][T];
        D_limit = Capacity;
        ans = init_sol2 + 1;
        
    }   else {
        T_limit = dmin[trans[S]][0][T];
        D_limit = tmin[trans[S]][1][T];
        ans = min(init_sol1, init_sol2) + 1;
    }
    
    // log_out << "PreSearch : " << ans << endl;
    // cout << tmin[trans[S]][T][1] << " " << dmin[trans[S]][T][1] << endl;
    
    LABEL init1{0, 0, S, 0};
    LABEL init2{0, 0, T, 1};
    priority_queue<LABEL, vector<LABEL>, greater<LABEL>> PQS;
    priority_queue<LABEL, vector<LABEL>, greater<LABEL>> PQT;
    PQS.push(init1);
    PQT.push(init2);
    PushLabel(init1);
    PushLabel(init2);
    int PQSSIZE = 1;
    int PQTSIZE = 1;
    ID[0] = trans[S];
    ID[1] = trans[T];
    PQNUM = SNUM = 0;
    // memset(Tmin, -1, sizeof(Tmin));
    // memset(Lmin, -1, sizeof(Lmin));
    fill(Tmin[0].begin(), Tmin[0].end(), -1);
    fill(Tmin[1].begin(), Tmin[1].end(), -1);
    fill(Lmin[0].begin(), Lmin[0].end(), -1);
    fill(Lmin[1].begin(), Lmin[1].end(), -1);
    
    while (PQSSIZE || PQTSIZE) {
        LABEL ulabel;
        int flag;
        if ((PQSSIZE < PQTSIZE && PQSSIZE != 0) || PQTSIZE == 0) {
            ulabel = PQS.top(); PQS.pop(); PQSSIZE--; flag = 0;
        }   else {
            ulabel = PQT.top(); PQT.pop(); PQTSIZE--; flag = 1;
        }
        
        bool test = true;
        for (auto &that_label : skyline[ulabel.index][ulabel.pos]) {
            if (ulabel == that_label) {
                test = false;
                break;
            }
        }
        if (test) continue;
        
        // AAAI的优化
        if (ulabel.sum_t >= Tmin[ulabel.index][ulabel.pos] && Tmin[ulabel.index][ulabel.pos] != -1) continue;
        Tmin[ulabel.index][ulabel.pos] = ulabel.sum_t;
        
        if (LabelPruning2(ulabel, ans)) continue;
        
        // 区间优化
        int newL = F(ulabel.sum_t, ulabel.sum_d + dmin[ID[1 - ulabel.index]][1][ulabel.pos]);
        if (newL >= Lmin[ulabel.index][ulabel.pos] && Lmin[ulabel.index][ulabel.pos] != -1) continue;
        Lmin[ulabel.index][ulabel.pos] = newL;
        
        /* CSP Pruning
        if (ulabel.sum_t + tmin[ID[1 - ulabel.index]][ulabel.pos][0] < ans &&
            ulabel.sum_d + tmin[ID[1 - ulabel.index]][ulabel.pos][1] <= R_limit)
            ans = ulabel.sum_t + tmin[ID[1 - ulabel.index]][ulabel.pos][0];
        if (ulabel.sum_t + dmin[ID[1 - ulabel.index]][ulabel.pos][0] < ans &&
            ulabel.sum_d + dmin[ID[1 - ulabel.index]][ulabel.pos][1] <= R_limit)
            ans = ulabel.sum_t + dmin[ID[1 - ulabel.index]][ulabel.pos][0];
        */
        
        if (ulabel.sum_d <= D_limit / 2) {
            LABEL vlabel;
            for (auto e : graph[ulabel.pos]) {
                vlabel.pos = e.v;
                vlabel.sum_t = ulabel.sum_t + e.t;
                vlabel.sum_d = ulabel.sum_d + e.d;
                vlabel.index = ulabel.index;
                
                if (LabelPruning2(vlabel, ans)) continue;
                
                if (!PushLabel(vlabel)) continue;
                
                if (flag == 0) {
                    PQS.push(vlabel); PQSSIZE++;
                }   else {
                    PQT.push(vlabel); PQTSIZE++;
                }
                // pre[vlabel] = ulabel;

                //pareto_optimal2[vlabel.index][vlabel.pos].Insert(vlabel.sum_t, vlabel.sum_d);
                PQNUM++;
            }
            // skyline[ulabel.index][ulabel.pos].push_back(ulabel);
            // dirt.push(make_pair(ulabel.index, ulabel.pos));
        }
        
        if (dmin[ID[!ulabel.index]][1][ulabel.pos] > D_limit / 2) continue;
        
        LabelCoupling2(ulabel, ans);
    }
    
    while (!dirt.empty()) {
        auto pr = dirt.front(); dirt.pop();
        skyline[pr.first][pr.second].clear();
    }
}

// 第二步 即从中继节点u到T，用Dijkstra解决
// int second_dis[MAXN], second_dis_t[MAXN], second_dis_d[MAXN];
// EDGE dijk_pre[MAXN];
vector<int> second_dis(MAXN), second_dis_t(MAXN), second_dis_d(MAXN);
vector<EDGE> dijk_pre(MAXN);
void SecondStep() {
    // log_out << "--------   Second Step Part   -------" << endl;
    
    // memset(second_dis, 0x3f, sizeof(second_dis));
    // memset(second_dis_t, 0, sizeof(second_dis_t));
    // memset(second_dis_d, 0, sizeof(second_dis_d));
    fill(second_dis.begin(), second_dis.end(), INF);
    fill(second_dis_t.begin(), second_dis_t.end(), 0);
    fill(second_dis_d.begin(), second_dis_d.end(), 0);
    second_dis[TT] = 0;
    
    priority_queue<pair<int,int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    pq.push(make_pair(second_dis[TT], TT));
    
    if (TT_is_station == false) {
        for (auto st : stations) {
            if (st == SS && SS_is_station == false) continue;
            if (st == TT) continue;
            Solve2(st, TT);
            if (ans_t == -1 || ans_d == -1) continue;
            second_dis[st] = ans;
            second_dis_t[st] = ans_t;
            second_dis_d[st] = ans_d;
            dijk_pre[st] = EDGE{TT, st, ans_t, ans_d};
            pq.push(make_pair(second_dis[st], st));
            log_out << st << " => " << TT;
            log_out << " ans=(" << ans_t << "," << ans_d << ")";
            log_out << " tmin=(" << tmin[trans[st]][0][TT] << "," << tmin[trans[st]][1][TT] << ")";
            log_out << " dmin=(" << dmin[trans[st]][0][TT] << "," << dmin[trans[st]][1][TT] << ")" << endl;
        }
    }
    
    while (!pq.empty()) {
        pair<int, int> pr = pq.top(); pq.pop();
        if (pr.first != second_dis[pr.second]) continue;
        
        // log_out << pr.second << endl;
        
        for (auto e : sketch[pr.second]) {
            int w = F(e.t, e.d);
            if (second_dis[e.u] + w < second_dis[e.v]) {
                second_dis[e.v] = second_dis[e.u] + w;
                pq.push(make_pair(second_dis[e.v], e.v));
                
                second_dis_t[e.v] = second_dis_t[e.u] + e.t;
                second_dis_d[e.v] = second_dis_d[e.u] + e.d;
                dijk_pre[e.v] = e;
            }
        }
    }
}



int D_upperbound[2], mid_point;
int ans_t1, ans_t2, ans_d1, ans_d2;
/*
从左往右
t_max ---------------- t_min
d_min ---------------- d_max
*/


void LabelCoupling1(LABEL &label, int &ans, int T, int after) {
    if (skyline[!label.index][label.pos].size()) {
        // int tmin = skyline[!label.index][label.pos].back().sum_d;
        // if (label.sum_t + tmin > T_limit) return true;
        
        for (auto that_label : skyline[1 - label.index][label.pos]) {
            SNUM++;
            // T_limit和D_limit无效
            // if (label.sum_t + that_label.sum_t > T_limit) continue;
            // if (label.sum_d + that_label.sum_d > D_limit) break;
            
            if (label.sum_d + that_label.sum_d < initial_d) continue;
            // if (ans <= F(label.sum_t + tmin, label.sum_d + that_label.sum_d - initial_d)) break;
            // if (total_ans <= F(label.sum_t + tmin, label.sum_d + that_label.sum_d - initial_d) + after) break;

            int subans = F(label.sum_t + that_label.sum_t, label.sum_d + that_label.sum_d - initial_d);
            if (subans + after < total_ans) {
                ans = subans;
                total_ans = subans + after;
                ans_label1 = label;
                ans_label2 = that_label;
                ans_t1 = label.sum_t;
                ans_t2 = that_label.sum_t;
                ans_d1 = label.sum_d;
                ans_d2 = that_label.sum_d;
                if (label.index == 1) { swap(ans_t1, ans_t2); swap(ans_d1, ans_d2); }
                mid_point = label.pos;
                // log_out << FindIndex(label) << " " << FindIndex(that_label) << endl;
            }
        }
    }
}

bool LabelPruning1(LABEL &label, int &ans, int after) {
    int dmin_t, dmin_d, tmin_t, tmin_d;
    
    // 到T的最省电的路线都已经超过Capacity
    if (dmin[ID[!label.index]][1].count(label.pos) == 0) {
        return true;
    }   else {
        dmin_t = dmin[ID[!label.index]][0][label.pos];
        dmin_d = dmin[ID[!label.index]][1][label.pos];
    }
    
    // 如果到T的最短路线超过了Capacity，则无法使用这条道路进行pruning
    if (tmin[ID[!label.index]][1].count(label.pos) == 0) {
        tmin_t = tmin_d = 0;
    }   else {
        tmin_t = tmin[ID[!label.index]][0][label.pos];
        tmin_d = tmin[ID[!label.index]][1][label.pos];
    }
    
    // 到终点的距离大于两个上限的和
    if (label.sum_d + dmin_d > D_upperbound[0] + D_upperbound[1]) return true;
    
    // 单边的上限
    if (label.sum_d > D_upperbound[label.index]) return true;
    
    // 预估的F值加上后继路径已经超过了当前的最优解
    if (F(label.sum_t + tmin_t,
          label.sum_d + dmin_d - initial_d)
        + after >= total_ans) return true;
    
    // 当前的t/d已经大于贪心上界
    // 属于支配性判定 可以继续使用
    if (dmin[ID[label.index]][0].count(label.pos) == 1) {
        if (label.sum_t > dmin[ID[label.index]][0][label.pos]) return true;
    }
    if (tmin[ID[label.index]][1].count(label.pos) == 1) {
        if (label.sum_d > tmin[ID[label.index]][1][label.pos]) return true;
    }
    
    // 同样属于支配性判定 可以继续使用
    if (label.sum_t > Tmin[label.index][label.pos] && Tmin[label.index][label.pos] != -1) return true;
    
    int newL = F(label.sum_t, label.sum_d + dmin[ID[!label.index]][1][label.pos] - initial_d);
    if (newL > Lmin[label.index][label.pos] && Lmin[label.index][label.pos] != -1) return true;
    
    return false;
}

void Solve1(int S, int T, int after) {
    if (S == T) { ans = F(0, 0); ans_t = ans_d = 0; return; }
    
    ans = INT_MAX;
    ans_t = ans_d = -1;
    
    LABEL init1 = {0, 0, S, 0};
    LABEL init2 = {0, 0, T, 1};
    PushLabel(init1);
    PushLabel(init2);
    priority_queue<LABEL, vector<LABEL>, greater<LABEL>> PQS;
    priority_queue<LABEL, vector<LABEL>, greater<LABEL>> PQT;
    PQS.push(init1);
    PQT.push(init2);
    int PQSSIZE = 1;
    int PQTSIZE = 1;
    ID[0] = trans[S];
    ID[1] = trans[T];
    PQNUM = SNUM = 0;
    // memset(Tmin, -1, sizeof(Tmin));
    // memset(Lmin, -1, sizeof(Lmin));
    fill(Tmin[0].begin(), Tmin[0].end(), -1);
    fill(Tmin[1].begin(), Tmin[1].end(), -1);
    fill(Lmin[0].begin(), Lmin[0].end(), -1);
    fill(Lmin[1].begin(), Lmin[1].end(), -1);
    
    // cout << "---------------------" << endl;
    
    while (PQSSIZE || PQTSIZE) {
        LABEL ulabel;
        int flag;
        if ((PQSSIZE < PQTSIZE && PQSSIZE != 0) || PQTSIZE == 0) {
            ulabel = PQS.top(); PQS.pop(); PQSSIZE--; flag = 0;
        }   else {
            ulabel = PQT.top(); PQT.pop(); PQTSIZE--; flag = 1;
        }

        bool test = true;
        for (auto &that_label : skyline[ulabel.index][ulabel.pos]) {
            if (ulabel == that_label) {
                test = false;
                break;
            }
        }
        if (test) continue;
        
        // cout << ulabel.index << " " << ulabel.pos << " " << ulabel.sum_t << " " << ulabel.sum_d << endl;
        
        // AAAI的优化
        // dominate的判定是必然成立的 如果被支配那必然可以被去掉
        if (ulabel.sum_t >= Tmin[ulabel.index][ulabel.pos] && Tmin[ulabel.index][ulabel.pos] != -1) continue;
        Tmin[ulabel.index][ulabel.pos] = ulabel.sum_t;
        
        if (LabelPruning1(ulabel, ans, after)) continue;
        
        // 区间优化
        int newL = F(ulabel.sum_t, ulabel.sum_d + dmin[ID[!ulabel.index]][1][ulabel.pos] - initial_d);
        if (newL > Lmin[ulabel.index][ulabel.pos] && Lmin[ulabel.index][ulabel.pos] != -1) continue;
        Lmin[ulabel.index][ulabel.pos] = newL;
        
        /* CSP Pruning
        if (ulabel.sum_t + tmin[ID[1 - ulabel.index]][ulabel.pos][0] < ans &&
            ulabel.sum_d + tmin[ID[1 - ulabel.index]][ulabel.pos][1] <= R_limit)
            ans = ulabel.sum_t + tmin[ID[1 - ulabel.index]][ulabel.pos][0];
        if (ulabel.sum_t + dmin[ID[1 - ulabel.index]][ulabel.pos][0] < ans &&
            ulabel.sum_d + dmin[ID[1 - ulabel.index]][ulabel.pos][1] <= R_limit)
            ans = ulabel.sum_t + dmin[ID[1 - ulabel.index]][ulabel.pos][0];
        */

        // 1/2的约束改为各自方向的约束
        if (ulabel.sum_d <= D_upperbound[ulabel.index]) {
            LABEL vlabel;
            for (auto e : graph[ulabel.pos]) {
                vlabel.pos = e.v;
                vlabel.sum_t = ulabel.sum_t + e.t;
                vlabel.sum_d = ulabel.sum_d + e.d;
                vlabel.index = ulabel.index;
                
                // cout << vlabel.index << " " << vlabel.pos << " " << vlabel.sum_t << " " << vlabel.sum_d << endl;
                
                if (LabelPruning1(vlabel, ans, after)) continue;
                
                if (!PushLabel(vlabel)) continue;
                
                if (flag == 0) {
                    PQS.push(vlabel); PQSSIZE++;
                }   else {
                    PQT.push(vlabel); PQTSIZE++;
                }

                PQNUM++;
            }
            // CheckDominance(ulabel);
            // skyline[ulabel.index][ulabel.pos].push_back(ulabel);
            // dirt.push(make_pair(ulabel.index, ulabel.pos));
        }
        
        // cout << "--------------------------" << endl;
        
        // 1/2的约束改为各自方向的约束
        if (dmin[ID[!ulabel.index]][1][ulabel.pos] > D_upperbound[!ulabel.index]) continue;
        
        // 在充电桩处才能进行解的组合
        if (is_station[ulabel.pos] == 1) LabelCoupling1(ulabel, ans, T, after);
    }
    
    while (!dirt.empty()) {
        auto pr = dirt.front(); dirt.pop();
        skyline[pr.first][pr.second].clear();
        //pareto_optimal2[pr.first][pr.second].Clear();
    }
}


// 第一步，即从S出发到中继节点u，u需要枚举
// too far表示S无法到u
// too close表示S出发到u后电量还是太多（因为需要从S到u，充电后继续往后走）
void FirstStep(int S) {
    // log_out << "--------   First Step Part   -------" << endl;
    
    total_ans = INT_MAX;
    
    for (int i = 0;i < stations.size(); i++) {
        int T = stations[i];
        int after = second_dis[stations[i]];
        if (after == 0x3f3f3f3f) continue;
        
        D_upperbound[0] = initial_d;
        D_upperbound[1] = Capacity;
        mid_point = -1;
        
        Solve1(S, T, after);
        
        // log_out << "PQNUM = " << PQNUM << endl;
        // log_out << "SNUM = " << SNUM << endl;
        
        if (mid_point != -1) {
            // log_out << S << " => " << mid_point << " => " << T << " (" << ans_t1 << "," << ans_d1 << ") (" << ans_t2 << "," << ans_d2 << ") F=" << ans << endl;
            // log_out << T;
            // while (T != TT) { T = dijk_pre[T].u; log_out << " => " << T; }
            // log_out << " F=" << after << endl;
            // log_out << "Total F = " << ans + after << endl;
            total_ans = min(total_ans, ans + after);
            // log_out << "----------------" << endl;
        }   else {
            // log_out << S << " => " << T << " =>* " << TT << " skip " << endl;
            // log_out << "----------------" << endl;
        }
    }
    log_out << "Final ans = " << total_ans << endl;
}

bool check(int u, int d) {
    for (auto &v : stations) {
        if (v == SS || v == TT) continue;
        if (dmin[trans[v]][1][u] <= d) return true;
    }
    return false;
}

int main() {
    // srand(6666);
    Init();
    ReadGraph();
    PreProcess(0, station_num - 1);
    
    SS_is_station = TT_is_station = false;
    for (int enum_d = 0.1 * Capacity; enum_d <= 0.5 * Capacity; enum_d += 0.1 * Capacity) {
        initial_d = enum_d;
        
        for (int enum_ST = 0; enum_ST < 2; enum_ST++) {

            SS = rand() % n + 1; while (is_station[SS] || !check(SS, initial_d)) SS = rand() % n + 1;
            TT = rand() % n + 1; while (is_station[TT] || !check(TT, Capacity)) TT = rand() % n + 1;
            //将SS和TT视作充电桩节点，然后进行导航
            trans[SS] = station_num + 1;
            trans[TT] = station_num + 2;
            stations[station_num] = SS;
            stations[station_num + 1] = TT;
            PreProcess(station_num, station_num + 1);

            log_out << "-------------------------" << endl;
            log_out << " SS = " << SS << " TT = " << TT << " initial_d = " << initial_d << endl;
            
            clock_t t1 = clock();
            
            SecondStep();
            
            clock_t t2 = clock();
            log_out << " SecondStep Time : " << (double)(t2 - t1) / CLOCKS_PER_SEC << " s  " << endl;
            
            FirstStep(SS);
            
            clock_t t3 = clock();
            log_out << " FirstStep Time : " << (double)(t3 - t2) / CLOCKS_PER_SEC << " s  " << endl;
            
            /*
            if (total_ans == INT_MAX) {
                log_out << "skip" << endl;
                enum_ST--;
                continue;
            }
             */
            
            exp_out << total_ans << " " << initial_d << " " << (double)(t2 - t1) / CLOCKS_PER_SEC << " " << (double)(t3 - t2) / CLOCKS_PER_SEC << endl;
        }
    }
}

/*
 1 10000
 ans = 2739556 (768716, 514168)
 
 1 100000 => Pre=4132926
 ans = 3777724 (1412755, 964234)
 4.6s PQ=1808745
 */

/*
16808 153722
 */
