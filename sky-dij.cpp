#include <iostream>
#include <vector>
#include <queue>
#include <time.h>
#include <limits.h>
#include <map>
#include <unordered_map>
#include <cmath>
#include <fstream>
#include <assert.h>
using namespace std;

typedef pair<int, int> pii;

const int MAXN = 2800000;
const int MAXM = 51 + 2;
const int INF = 0x3f3f3f3f;
const int Capacity = 2000000;

int n, m;
int station_num;

//自己是u,邻居是v,t是weight,d是cost
struct EDGE {
    int u, v, t, d;
};

vector<vector<EDGE>> graph(MAXN, vector<EDGE>());
vector<vector<EDGE>> sketch(MAXN, vector<EDGE>());//选取的50个充电桩
vector<vector<int>> label[MAXN];//store paths from start to vertex
vector<int> stations;
unordered_map<int, int> trans;//节点x对应第trans[x]个充电桩
unordered_map<int, int> is_station;
int is_visited[51];//表示某个节点到其他所有节点的skyline都算出来了
vector<vector<vector<int>>> station_path(51,vector<vector<int>>(51,vector<int>()));//station_path[i][j]表示充电桩i到j的最短路径

ifstream road_data("dataset/graph/USA-road-t-d.FLA.gr");
ifstream sketch_data("dataset/FLA-road-station.txt");

ofstream log_out("log/test-log-FLA.txt");
ofstream exp_out("log/test-exact-FLA.txt");

//充电函数，这里加上了travel time，也就是sum_t
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

void clear_label() {
    for (auto& vec : label) {
        vec.clear();
    }
}

void ReadGraph() {
    long long int num = 0;
    std::ios::sync_with_stdio(false);

    road_data >> n >> m;
    cout<<n<<" "<<m<<endl;
    for (int i = 1;i <= m; i++) {
        EDGE e;
        road_data >> e.u >> e.v >> e.t >> e.d;
//        if(num < 5 || num > m-5){
//            cout<< "num: "<< num << "  " << e.u<< e.v << e.t << e.d<< endl;
//        }
        graph[e.u].push_back(e);
        num++;
    }


    sketch_data >> station_num;
    for (int i = 1;i <= station_num; i++) {
        int x;
        sketch_data >> x;
        stations.push_back(x);
        trans[x] = i;
        is_station[x] = 1;
//        if(i <= 9) is_visited[i] = 1;
    }


    for (int i = 1; i <= station_num * (station_num - 1); i++) {
        EDGE e;
        sketch_data >> e.u >> e.v >> e.t >> e.d;
        if (e.t == -1 || e.d == -1) continue;
        sketch[e.u].push_back(e);
    }

    // 占位
//    stations.push_back(0);
//    stations.push_back(0);
    // task_data >> SS >> TT;
    // if (trans[SS] == 0) { stations.push_back(SS); trans[SS] = station_num + 1; SS_is_station = false; }
    // if (trans[TT] == 0) { stations.push_back(TT); trans[TT] = station_num + 2; TT_is_station = false; }
}

//todo:通过记录每个label的前缀，来找出一条完整的路径
void calculate_station_path(int start, int target){
//    for(auto target: stations){
//        if(station_path[start][target].empty()){
//        }
//    }
    int total_time = INT_MAX;
    vector<int> shortest_path;
    log_out<<1<<endl;
    for(auto p: label[target]){
        log_out<<2<<endl;
        int time = F(p[0], p[1]);
        if(time < total_time){
            total_time = time;
            shortest_path = p;
        }
    }
    station_path[start][target] = shortest_path;
}

void print(int index){
    for(auto p: label[index]){
        cout<<"sum_t: "<<p[0]<<" sum_d: "<<p[1]<<endl;
        cout<<"path: ";
        for(int i=2;i<p.size();i++){
            cout<<p[i]<<' ';
        }
        cout<<endl<<endl;
    }
}


bool is_dominating(vector<int> path1, vector<int> path2){
    if(path1[0]<=path2[0] && path1[1]<=path2[1]) return true;
    return false;
}

struct CompareVector {
    bool operator()(const vector<int>& v1, const vector<int>& v2) {
        return v1[1] > v2[1];  // 按照vector的第二个元素升序排列
    }
};


//todo: 如果遍历到充电桩,那么需要将electricity转换的时间先计算，然后从0重新开始计算,因为在充电桩可以重新充电,那么上一个充电桩到达这个充电桩的电量就为0(不需要做，因为在充电网络会考虑到这种情况)
//path[0]:sum_t,i.e.,time. path[1]:sum_d,i.e.,distance(electricity), path[2]: start(去掉了，因为起点肯定是s), path[3]:end
//NY:限制180万和250万次迭代
//COL:限制30万次
//FLA:限制350万次
bool sky_dijkstra(int s,int t){
    int num = 0;
    priority_queue<vector<int>, vector<vector<int>>, CompareVector> q;
    q.push({0,0,s});
    label[s].push_back({0,0,s});
    while(!q.empty()){
        vector<int> path=q.top(); q.pop();
        int path_end = path.back();
        if(path_end != t){
            for(auto edge : graph[path_end]){
                num++;
//                if(num>60) return;
                int expand = edge.v, flag = 1;
                vector<int> expand_path=path;
                expand_path.back() = expand;
                expand_path[0] += edge.t;
                expand_path[1] += edge.d;
                if(num > 3500000) return false;
//                if(num > 6000000){
//                    log_out<<"original path: ";
//                    for(auto v: path){
//                        log_out<<v<<' ';
//                    }
//                    log_out << " expand_vertex: "<<expand<<endl;
//                    log_out << "========="<<endl;
//                }
                if(expand_path[1] > Capacity) continue;
//                if(expand == 32715){
//                    log_out<<" segmentation fault: "<<endl;
//                    for(auto p1: label[expand]){
//                        for(auto v:p1){
//                            log_out<<v<<' ';
//                        }
//                        log_out<<endl;
//                    }
//                }
                if(label[expand].empty()){
//                    log_out<<" push expand path to pq"<<endl;
                    q.push(expand_path);
                    label[expand].push_back(expand_path);
                    continue;
                }
//                int ccc=0;
                for(auto p : label[expand]){
//                    log_out<<ccc++<<endl;
                    if(is_dominating(p, expand_path)){
//                        log_out<<"expand path is dominated"<<endl;
                        flag = 0; break;
                    }
                }
//                log_out<<" 1111"<<endl;
                if(flag == 1) {// 没有被dominated
                    vector<vector<int>> new_label;
                    for(auto p: label[expand]){//删除所有被expand_path dominate的path
                        if(is_dominating(expand_path, p)) continue;
                        new_label.push_back(p);
//                        log_out<<ccc++<<endl;
                    }
                    label[expand] = new_label;
//                    log_out<<"2222"<<endl;
//                    log_out<<" push expand path to pq"<<endl;
                    q.push(expand_path);
                    label[expand].push_back(expand_path);
                }
//                log_out << "===========" <<endl;
            }
        }else{
            return true;
        }
    }
}



int main(){
//    sky_dijkstra(1,5);
//    clock_t t1 = clock();
//    clock_t t2 = clock();
//    exp_out << "start: "<< 1 << "  target: " << 2 << "  distance: " << graph[1][2].d << "  construct time: " << (double)(t2 - t1) / CLOCKS_PER_SEC << endl;
//    exp_out << "start: "<< 1 << "  target: " << 2 << "  distance: " << graph[1][2].d;
    ReadGraph();
//    exp_out << "start: "<< 1 << "  target: " << 2 << "  distance: " << graph[1][2].d;
//    sky_dijkstra(1, 2);
//    calculate_station_path(s, t);
//    sky_dijkstra(101357,199443);

    for(auto s: stations){
        if(is_visited[trans[s]]) continue;
        for(auto t: stations){
            if(t == s || is_visited[trans[t]]) continue;
            clear_label();

            clock_t t1 = clock();

            bool solved = sky_dijkstra(s, t);
            log_out << "sky-dij finish" << endl;

//            calculate_station_path(1, 2);
//            cout <<" calculate_finish" <<endl;

            clock_t t2 = clock();

            int distance = 0;
            for(auto e: sketch[s]){
                if(e.v == t){
                    distance = e.d;
                }
            }
            if(solved){
                int min_d = INT_MAX;
                for(auto p: label[t]){
                    if(p[1] < min_d){
                        min_d = p[1];
                    }
                }
                distance = min_d;
            }

            double construct_time = solved ? (double)(t2 - t1) / CLOCKS_PER_SEC : (-1.0);
            exp_out << "start: "<< s << "  target: " << t << "  distance: " << distance << "  construct time: " << construct_time;
            exp_out << endl;
        }
        is_visited[trans[s]] = 1;
    }
}