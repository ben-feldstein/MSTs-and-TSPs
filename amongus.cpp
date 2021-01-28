// Project Identifier: 9B734EC0C043C5A836EA0EBE4BEFEA164490B2C7
#include <getopt.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <iomanip>
#include <limits>

using namespace std;

class AmongUs {
public:
    void get_options(int argc, char *argv[]); //set mode
    void read_input(); //read in coordinates
    void mst_mode(); //build mst
    void fast_tsp_mode(); //build tsp estimate quickly
    void opt_tsp_mode(); //build optimal tsp
    void choose_algo(); //pick algorithm to use
    void output(); //output resulting path

private:
    struct Prim {
        double d = numeric_limits<double>::infinity(); //minimal edge weight to v
        int p = -1; //index of preceeding room
        bool k = false; //visited?
    };
    struct Salesman {
        double distance = 0;
        int prev = 0;
        int next = 0;
        bool visited = false;
    };
    enum class Location : uint8_t {
        Lab,
        Outer,
        Decont
    };
    struct Coordinate {
        Coordinate(int x_in, int y_in) {
            x = x_in;
            y = y_in;
        }
        int x;
        int y;
        Location zone;
    };
    enum class Mode : uint8_t {
        MST,
        FASTTSP,
        OPTTSP
    };
    vector<Coordinate> rooms; //input coordinates
    vector<Prim> mst; //minimum spanning tree container
    vector<Salesman> tsp; //fast tsp solution container
    vector<vector<double>> distance_matrix;
    vector<int> path;
    vector<int> best_path;
    double total_distance = numeric_limits<double>::infinity();
    double upper_bound = 0;
    double running_total = 0;
    int num_rooms = 0; //number of rooms
    Mode mode; //game mode
    double distance(const Coordinate& p1, const Coordinate& p2) { //calculate euclidean distance between two points
        if (mode == Mode::MST) {
            if ((p1.zone == Location::Lab && p2.zone == Location::Outer) ||
                (p2.zone == Location::Lab && p1.zone == Location::Outer)) {
                return numeric_limits<double>::infinity(); //if MST mode and the points are in different zones, return infinity
            }
        }
       /* if (mode == Mode::OPTTSP) {
            if (p1.x == p2.x && p1.y == p2.y) { return 0; }
        }*/
        double diff1 = static_cast<double>(p1.x - p2.x);
        double diff2 = static_cast<double>(p1.y - p2.y);
        return sqrt((diff1 * diff1) + (diff2 * diff2));
        //return hypot(diff1, diff2);
    } //distance()
    void set_zone(Coordinate& p, bool& L, bool& D, bool& O) { //set the zone for the point
        if ((p.x == 0 && p.y < 1) || (p.y == 0 && p.x < 1)) {
            p.zone = Location::Decont;
            D = true;
        }
        else if (p.x < 0 && p.y < 0) {
            p.zone = Location::Lab;
            L = true;
        }
        else {
            p.zone = Location::Outer;
            O = true;
        }
    } //set_zone()
    double calc_total_weight() {
        double total = 0;
        for (auto& it : tsp) {
            total += it.distance;
        }
        return total;
    } //calc_total_weight()
    void calculate_distance_matrix(vector<vector<double>>& v) {
        for (int i = 0; i < num_rooms; i++) {
            for (int j = 0; j < num_rooms; j++) {
                v[i][j] = distance(rooms[i], rooms[j]);
            }
        }
    }// calculate_distance_matrix()
    double low_bound_cost(vector<Prim>& vec, size_t permlength) {
        double res = 0;
        for (size_t i = permlength; i < path.size(); i++) {
            res += vec[path[i]].d;
        }
        return res;
    }// low_bound_cost()
    bool promising(size_t permlength) {
        /*if (permlength == 1 || permlength == 2) { return true; }*/
        vector<Prim> low_bound(num_rooms);
        for (size_t i = 0; i < permlength; i++) {
            //low_bound[path[i]].d = 0;
            low_bound[path[i]].k = true;
        }
        low_bound[path[permlength]].d = 0;
        low_bound[path[permlength]].p = -1;
        int current = -1;
        int j = static_cast<int>(permlength);
        while (j <= num_rooms) {
            double min_distance = numeric_limits<double>::infinity();
            for (int i = static_cast<int>(permlength); i < num_rooms; i++) {
                if (low_bound[path[i]].k == false) {
                    if (low_bound[path[i]].d < min_distance) {
                        min_distance = low_bound[path[i]].d;
                        current = i;
                    }
                }
            } //for i
            low_bound[path[current]].k = true;
            for (int i = static_cast<int>(permlength); i < num_rooms; i++) {
                if (low_bound[path[i]].k == false) {
                    if (distance_matrix[path[current]][path[i]] < low_bound[path[i]].d) {
                        low_bound[path[i]].d = distance_matrix[path[current]][path[i]];
                        low_bound[path[i]].p = path[current];
                    }
                }
            } //for i
            j++;
        } //while j
        double front_arm = numeric_limits<double>::infinity();
        for (int i = static_cast<int>(permlength); i < num_rooms; i++) {
            if (distance_matrix[path[0]][path[i]] < front_arm) {
                front_arm = distance_matrix[path[0]][path[i]];
            }
        }
        double back_arm = numeric_limits<double>::infinity();
        for (int i = static_cast<int>(permlength); i < num_rooms; i++) {
            if (distance_matrix[path[permlength-1]][path[i]] < back_arm) {
                back_arm = distance_matrix[path[permlength-1]][path[i]];
            }
        }
        //cout << "PermLength: "<< permlength << " Running Total: " << running_total << " MST Est.: " << low_bound_cost(low_bound, permlength)<< " front arm: " << front_arm
        //     << " back arm: " << back_arm << " upper b: " << upper_bound << "\n"; //test
        if (front_arm + back_arm + low_bound_cost(low_bound, permlength) + running_total > upper_bound) {
            return false;
        }
        else {
            return true;
        }
    }// promising()
    void genPerms(size_t permLength) {
        if (permLength == path.size()) {
            double new_total = running_total + distance_matrix[0][path[permLength - 1]];
            if (new_total < upper_bound) {
                best_path = path;
                upper_bound = new_total;
            }
            return;
        } // if
        if (!promising(permLength)) {
            return;
        }
        for (size_t i = permLength; i < path.size(); ++i) {
            swap(path[permLength], path[i]);
            running_total += distance_matrix[path[permLength -1]][path[permLength]];
            genPerms(permLength + 1);
            running_total -= distance_matrix[path[permLength -1]][path[permLength]];
            swap(path[permLength], path[i]);
            
        } // for
    } // genPerms()
}; //AmongUs

int main(int argc, char *argv[]) {
    ios_base::sync_with_stdio(false);
    cout << std::setprecision(2);
    cout << std::fixed;
    //std::cout << std::setw(10); //test

    AmongUs game;
    game.get_options(argc, argv);
    game.read_input();
    game.choose_algo();
    game.output();

    return 0;
}// main()

void AmongUs::get_options(int argc, char *argv[]) {
    string mode_in = "";
    int option_index = 0, option = 0;
    opterr = true;

    struct option longOpts[] = { { "mode", required_argument, nullptr, 'm' },
                                { "help", no_argument, nullptr, 'h' },
                                { nullptr, 0, nullptr, '\0' } };
    if (argc == 1) {
        cerr << "Error: No mode specified\n";
        exit(1);
    }
    while ((option = getopt_long(argc, argv, "m:h", longOpts, &option_index)) != -1) {
        switch (option) {
        case 'm':
            mode_in = optarg;
            if (!optarg) {
                cerr << "Error: No mode specified\n";
                exit(1);
            }
            if (mode_in == "MST") {
                mode = Mode::MST;
            }
            else if (mode_in == "FASTTSP") {
                mode = Mode::FASTTSP;
            }
            else if (mode_in == "OPTTSP") {
                mode = Mode::OPTTSP;
            }
            else {
                cerr << "Error: Invalid mode\n";
                exit(1);
            }
            break;
        case 'h':
            std::cout << "This program takes command line arguments and a text file of coordinates,\n"
                << "it reads the coordinates and then completes a graph algorithm specified on the command line,\n"
                << "the program can create a minimum spanning tree, a fast traveling salesperson solution,\n"
                << "and an optimal traveling salesperson solution.\n"
                << "Usage: \'./amongus\n\t[--help | -h]\n"
                << "\t[--mode | -m < MST | FASTTSP | OPTTSP >]\n"
                << "\t< <coordinates .txt file>\'" << std::endl;
            exit(0);
        default:
            cerr << "Error: Invalid command line option" << endl;
            exit(1);
        } //switch
    } //while
} //get_options()

void AmongUs::read_input() {
    bool has_lab = false;
    bool has_decont = false;
    bool has_outer = false;
    cin >> num_rooms;
    rooms.reserve(num_rooms);
    for (int i = 0; i < num_rooms; i++) {
        int x, y;
        cin >> x >> y;
        Coordinate c(x, y);
        if (mode == Mode::MST) {
            set_zone(c, has_lab, has_decont, has_outer);
        }
        rooms.push_back(c);
    }
    if (mode == Mode::MST) {
        if (has_lab && has_outer && !has_decont) {
            cerr << "Cannot construct MST\n";
            //cout << has_lab << " " << has_outer << " " << has_decont << "\n";
            exit(1);
        }
    }
} //read_input()

void AmongUs::mst_mode() {
    mst.resize(num_rooms);
    mst[0].d = 0;
    mst[0].p = -1;
    int current = -1;
    int j = 0;
    while (j < num_rooms) {
        double min_distance = numeric_limits<double>::infinity();
        for (int i = 0; i < num_rooms; i++) {
            if (mst[i].k == false) {
                if (mst[i].d < min_distance) {
                    min_distance = mst[i].d;
                    current = i;
                }
            }
        } //for i
        mst[current].k = true;
        for (int i = 0; i < num_rooms; i++) {
            if (mst[i].k == false) {
                if (distance(rooms[current], rooms[i]) < mst[i].d) {
                    mst[i].d = distance(rooms[current], rooms[i]);
                    mst[i].p = current;
                }
            }
        } //for i
        j++;
    } //while j
} //mst_mode()

void AmongUs::fast_tsp_mode() {
    int num_visited = 0;
    tsp.resize(num_rooms);
    tsp[0].visited = true;
    num_visited++;
    double min_distance = numeric_limits<double>::infinity();
    int next = 0;
    int next_next = 0;
    double min_dis_ii = numeric_limits<double>::infinity();
    for (int i = 1; i < num_rooms; i++) {
        double dis = distance(rooms[0], rooms[i]);
        if (dis < min_distance) {
            min_distance = dis;
            next = i;
        }
    }// for i
    for (int i = 1; i < num_rooms; i++) {
        if (i != next) {
            double dis = distance(rooms[next], rooms[i]);
            if (dis < min_dis_ii) {
                min_dis_ii = dis;
                next_next = i;
            }
        } 
    }// for i
    tsp[0].next = next;
    tsp[0].prev = next_next;
    tsp[next].visited = true;
    num_visited++;
    tsp[next].next = next_next;
    tsp[next].prev = 0;
    tsp[next].distance = min_distance;
    total_distance = min_distance;
    tsp[next_next].visited = true;
    num_visited++;
    tsp[next_next].next = 0;
    tsp[next_next].prev = next;
    tsp[next_next].distance = min_dis_ii;
    total_distance += min_dis_ii;
    int k = 0;
    while (num_visited < num_rooms) {
        for (int i = k+1; i < num_rooms; i++) {
            if (tsp[i].visited == false) {
                k = i;
                break;
            }
        }// for i
        double min = numeric_limits<double>::infinity();
        int vertex_i = 0;
        for (int i = 0; i < num_rooms; i++) {
            if (tsp[i].visited) {
                double new_total = total_distance - distance(rooms[i], rooms[tsp[i].next]) +
                    distance(rooms[i], rooms[k]) + distance(rooms[k], rooms[tsp[i].next]);
                if (new_total < min) {
                    vertex_i = i;
                    min = new_total;
                }
            }
        }// for i
        int vertex_j = tsp[vertex_i].next;
        tsp[vertex_i].next = k;
        tsp[k].prev = vertex_i;
        tsp[k].next = vertex_j;
        tsp[k].visited = true;
        num_visited++;
        tsp[vertex_j].prev = k;
        tsp[k].distance = distance(rooms[vertex_i], rooms[k]);
        tsp[vertex_j].distance = distance(rooms[vertex_j], rooms[k]);
        total_distance = min;
    } //while
    tsp[0].distance = distance(rooms[0], rooms[tsp[0].prev]);
    total_distance += tsp[0].distance;
    if (mode == Mode::OPTTSP) {
        path.resize(num_rooms);
        path[0] = 0;
        int index = tsp[0].next;
        for (int i = 1; i < num_rooms; i++) {
            path[i] = index;
            index = tsp[index].next;
        }
        best_path = path;
    }
} //fast_tsp_mode()

void AmongUs::opt_tsp_mode() {
    fast_tsp_mode();
    upper_bound = calc_total_weight();
    distance_matrix.resize(num_rooms, vector<double>(num_rooms));
    calculate_distance_matrix(distance_matrix);
    genPerms(1);
} //opt_tsp_mode()

void AmongUs::choose_algo() {
    switch (mode) {
    case Mode::MST:
        mst_mode();
        break;
    case Mode::FASTTSP:
        fast_tsp_mode();
        break;
    case Mode::OPTTSP:
        opt_tsp_mode();
        break;
    } //switch
} //choose_algo()

void AmongUs::output() {
    if (mode == Mode::MST) {
        double total = 0;
        for (auto& it : mst) {
            total += it.d;
        }
        cout << total << "\n";
        for (int i = 1; i < num_rooms; i++) {
            cout << min(mst[i].p, i) << " " << max(mst[i].p, i) << " \n";
        }
    } //MST mode
    else if (mode == Mode::FASTTSP){
        cout << calc_total_weight() << "\n";
        cout << "0 ";
        int index = tsp[0].next;
        for (int i = 1; i < num_rooms; i++) {
            cout << index << " ";
            index = tsp[index].next;
        }
    } //FASTTSP mode
    else {
        cout << upper_bound << "\n";
        for (int i = 0; i < num_rooms; i++) {
            cout << best_path[i] << " ";
        }
    }// OPTTSP mode
} //output()