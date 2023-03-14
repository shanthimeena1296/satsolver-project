// Compile with c++ ece650-a2cpp -std=c++11 -o ece650-a2
#include <iostream>
#include <sstream>
#include <fstream>
#include <regex>
#include <vector>
#include <queue>
#include <cstdio>
#include <string>
#include <algorithm>
#include <memory>
#include "minisat/core/SolverTypes.h"
#include "minisat/core/Solver.h"
#include <pthread.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <chrono>

using namespace Minisat;

//declaring
std::vector < int > solveVertexCover(int, std::vector < int > , int);
std::vector < int > find_vertexcover(int, std::vector < int > [], int, int);

std::stringstream outVC, outVC1, outVC2;
std::stringstream outVC_time, outVC1_time, outVC2_time;

int v;
int edge_count;
bool g_cancel = 0;
int g_lock = 0;
bool g_stopwatchdog = 1;
bool g_RunOnce = 0;

void acquire(int * lock) {
  while ( * lock != 0);
  * lock = 1;
}

void release(int * lock) {
  * lock = 0;
}

// https://www.geeksforgeeks.org/shortest-path-unweighted-graph/

// Function to class graph for vert and edge
class create_graph {
  int vert;
  public:
    std::vector < int > * adj_l;
  explicit create_graph(int vert) {
    this -> vert = vert;
    this -> adj_l = new std::vector < int > [vert + 1];
  };

};
create_graph * undirected_graph;
// function to add edges
void edge_add(std::vector < int > vert_array[], int v1, int v2) {
  vert_array[v1].push_back(v2);
  vert_array[v2].push_back(v1);
}

void * watchdog(void * arg) {

  pthread_t t = * (pthread_t * ) arg;
  // std::cout<<t<<std::endl<<std::flush;

  while (!g_cancel) {
    if (g_stopwatchdog == 0) {
      std::time_t start = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
      std::time_t k = start + 2;
      while (start != k) {
        start = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        if (g_stopwatchdog == 1) {
          // std::cout<<"timer stopped before ending"<<std::endl;
          g_cancel = 1;
          break;
        }
      }
      if (start == k) {
        outVC << "CNF-SAT-VC: timeout" << std::endl;
        pthread_cancel(t);
        g_stopwatchdog = 1;
        g_cancel = true;

      }

    }

    // std::cout<<"5 seconds past";
  }
  // std::cout<<"watchdog ending..."<<std::endl;
  return nullptr;
}
//For vertex cover solving
int solveVertexCover(int vertex, std::vector < int > edge[], int edge_count) {

  std::unique_ptr < Minisat::Solver > solver(new Minisat::Solver());

  std::vector < int > value;
  std::vector < int > outValue;

  int vect_size = vertex;
  int min = 1;
  std::vector < int > tmp = {
    -1
  };
  //calculating the mid value
  while (vect_size >= min) {

    int mid_value = (vect_size + min) / 2;

    value = find_vertexcover(vertex, edge, mid_value, edge_count);
    //checking whether it sat or not
    bool nosat = std::equal(value.begin(), value.end(), tmp.begin());
    //If not satisfy, then do below

    if (not nosat) {
      vect_size = mid_value - 1;
      outValue.clear();
      outValue = value;
    } else {
      min = mid_value + 1;
    }
  }
  //sort the vector cover values in ascending
  sort(outValue.begin(), outValue.end());

  //Printing the final output here
  outVC << "CNF-SAT-VC: ";
  for (int i = 0; i < outValue.size(); ++i) {
    if (i == outValue.size() - 1) { //print like this for the last item
      outVC << outValue[i] << std::endl;
    } else {
      outVC << outValue[i] << ",";
    }

  }
  // std::cout<<out.str()<<std::flush;
  // std::cout << "\n";
  return outValue.size();
}

void * solveVertexCoverThread(void * arg) {

  clock_t ON, OFF;
  double duration;

  int test_v = * (int * ) arg;
  // std::cout<<pthread_self()<<std::endl<<std::flush;
  pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL); //allow async termination of thread

  // std::vector<int>* k; //USE THIS TO FIX STD::VECTOR<INT>* ERROR
  // while (!g_cancel){
  // if (g_RunOnce){
  ON = clock();
  // solveVertexCover(std::ref(v), std::ref(undirected_graph -> adj_l), std::ref(edge_count));
  // alarm(5);
  g_stopwatchdog = 0;
  solveVertexCover(test_v, undirected_graph -> adj_l, edge_count);
  g_stopwatchdog = 1;

  OFF = clock();

  duration = (double)(OFF - ON) / CLOCKS_PER_SEC;
  outVC_time << "vc time: " << duration;
  // g_RunOnce=0;
  // g_cancel=1;
  // }
  // }
  // std::cout<<"cnf thread ending..."<<std::endl;
  std::cout << "solvevertextcover thread----------" << duration << std::endl;
  return nullptr;

}

//https://git.uwaterloo.ca/ece650-1221/minisat-example/-/blob/master/ece650-minisat.cpp
//using minisat for finding the vertex cover
std::vector < int > find_vertexcover(int vertex_size, std::vector < int > * edge, int mid, int edge_count) {

  std::unique_ptr < Minisat::Solver > solver(new Minisat::Solver());

  Minisat::Lit arr[vertex_size][mid];

  for (int i = 0; i < vertex_size; i++) {
    for (int j = 0; j < mid; j++) {
      arr[i][j] = Minisat::mkLit(solver -> newVar());
    }
  }
  // At least (exactly only one) one vertex is the ith vertex in the vertex cover, i in [1,k]
  for (int i = 0; i < mid; i++) {
    Minisat::vec < Minisat::Lit > clause;
    for (int j = 0; j < vertex_size; j++) {
      clause.push(arr[j][i]);
    }
    solver -> addClause(clause);
    clause.clear();

  }
  // No one vertex can appear twice in a vertex cover
  for (int m = 0; m < vertex_size; m++) {
    for (int p = 0; p < mid; p++) {
      for (int q = p + 1; q < mid; q++) {
        solver -> addClause(~arr[m][p], ~arr[m][q]);
      }
    }
  }
  // No more than one vertex appears in the mth position of the vertex cover.
  for (int mth_pos = 0; mth_pos < mid; mth_pos++) {
    for (int p = 0; p < vertex_size; p++) {
      for (int q = p + 1; q < vertex_size; q++) {
        solver -> addClause(~arr[p][mth_pos], ~arr[q][mth_pos]);
      }
    }
  }
  // Every edge is incident to at least one vertex in the vertex cover
  for (int i = 1; i <= vertex_size; i++) {
    Minisat::vec < Minisat::Lit > clause;
    for (int l = 0; l < edge[i].size(); l++) {
      for (int m = 0; m < mid; m++) {
        int src = i - 1;
        int dst = edge[i][l] - 1;

        clause.push(arr[src][m]);
        clause.push(arr[dst][m]);
      }
      solver -> addClause(clause);
    }
    clause.clear();
  }

  bool res = solver -> solve();

  std::vector < int > vertx_data;
  if (res) {
    for (int i = 0; i < vertex_size; i++) {
      for (int j = 0; j < mid; j++) {
        if (solver -> modelValue(arr[i][j]) == l_True) {
          vertx_data.push_back(i + 1);
        }
      }
    }
    return vertx_data;
  } else {
    return {
      -1
    };
  }
}

void print_vecvec(std::vector < std::vector < int >> vecvec) {
  std::stringstream out;

  for (int vert = 0; vert < vecvec.size(); vert++) {
    for (int k: vecvec[vert]) {
      out << "vert: " << vert << " | " << k << std::endl;
    }
  }
  std::cout << out.str() << std::flush;
}

std::vector < std::vector < int >> ArrayVec_To_VecVec(std::vector < int > arrayvec[], int vertex_size) {
  //bridge the work between MEENA and BRIAN work
  std::vector < std::vector < int >> vecvec;
  vecvec.resize(vertex_size + 1);

  for (int vertex = 0; vertex < vertex_size + 1; vertex++) { //converting the array of vecs to vecOfVec

    if (arrayvec[vertex].size() != 0) {
      for (int vert_neighbor: arrayvec[vertex]) {
        // std::cout<<"got here"<<vertex<<":"<<vert_neighbor<<std::endl;
        vecvec[vertex].push_back(vert_neighbor);
      }
    }
  }
  return vecvec;
}

int APPROX_VC_1(std::vector < int > edge[], int vertex_size) {
  // Taking the highest degree vertex and removing children
  int highest_degree = 0;
  int highest_degree_vert = 0;
  std::vector < int > assignments;
  bool sibling_info_empty = false;
  // std::stringstream outVC1;

  std::vector < std::vector < int >> siblings;

  siblings = ArrayVec_To_VecVec(edge, vertex_size);

  while (!sibling_info_empty) {
    // PrintSiblingDatabase(siblings);
    for (std::vector < int > vert: siblings) {
      if (vert.size() == 0) {
        sibling_info_empty = true;

      } else {
        sibling_info_empty = false;
        break;
      }
    }

    if (sibling_info_empty == true) {
      break;
    }

    for (int vert = 1; vert < siblings.size(); vert++) {
      if (siblings[vert].size() > 0) {
        if (siblings[vert].size() > highest_degree) {
          highest_degree = siblings[vert].size();
          highest_degree_vert = vert; //update the highest degree vert 
        }
      }
    }

    assignments.push_back(highest_degree_vert);
    // std::cout<<"Highest degree vert= "<<highest_degree_vert<<std::endl;
    for (int verts_to_erase: siblings[highest_degree_vert]) { //remove highest_degree vertex from all othe rvertex's info
      for (int i = 0; i < siblings[verts_to_erase].size(); i++) {
        if (siblings[verts_to_erase][i] == highest_degree_vert) {
          siblings[verts_to_erase].erase(siblings[verts_to_erase].begin() + i);
        }
      }
    }

    siblings[highest_degree_vert].clear();
    // print_vecvec(siblings);
    highest_degree = 0;

  }

  outVC1 << "APPROX-VC1: ";
  // for (int k : assignments){
  //     std::cout<<k<<",";
  // }

  std::sort(assignments.begin(), assignments.end()); //sort the sat assignments from lowest to highest

  for (int i = 0; i < assignments.size(); i++) { //print out the solution
    if (i == assignments.size() - 1) {
      outVC1 << assignments[i] << std::endl;
    } else {
      outVC1 << assignments[i] << ",";
    }

  }
  // std::cout<<out.str()<<std::flush;
  return assignments.size();
}
//approx vc thread time calculation
void * APPROX_VC_1Thread(void * arg) {

  clock_t ON, OFF;
  double duration;
  int test_v = * (int * ) arg;
  // std::cout<<pthread_self()<<std::endl<<std::flush;
  // while (!g_cancel){
  //     if (g_RunOnce){
  ON = clock();
  APPROX_VC_1(std::ref(undirected_graph -> adj_l), test_v);

  OFF = clock();

  duration = (double)(OFF - ON) / CLOCKS_PER_SEC;
  outVC1_time << "vc1 time: " << duration;
  //     }
  // }
  std::cout << "approxVCThread 1----------" << duration << std::endl;
  return nullptr;

}

int APPROX_VC_2(std::vector < int > edge[], int vertex_size) {
  // Taking the highest degree vertex and removing children
//   int highest_degree = 0;
//   int highest_degree_vert = 0;
  std::vector < int > assignments;
  bool sibling_info_empty = false;

  // std::stringstream outVC2;
  int u, v;

  std::vector < std::vector < int >> siblings;

  siblings = ArrayVec_To_VecVec(edge, vertex_size);

  while (!sibling_info_empty) {
    // PrintSiblingDatabase(siblings);
    for (std::vector < int > vert: siblings) {
      if (vert.size() == 0) {
        sibling_info_empty = true; //checking if all siblings are covered

      } else {
        sibling_info_empty = false;
        break;
      }
    }

    if (sibling_info_empty == true) {
      break;
    }
    // print_vecvec(siblings);
    // std::cout<<"------"<<std::endl;
    //pick an edge to use

    for (int vert = 1; vert < siblings.size(); vert++) {
      if (siblings[vert].size() > 0) {
        u = vert;

        v = siblings[vert].back();
        siblings[vert].pop_back();
        break;
      }
    }

    assignments.push_back(v);
    assignments.push_back(u);
    // std::cout<<u<<" : "<<v<<std::endl;

    // for(int verts_to_erase :siblings[highest_degree_vert]){ //remove highest_degree vertex from all othe rvertex's info
    //     for(int i=0; i<siblings[verts_to_erase].size(); i++){
    //         if (siblings[verts_to_erase][i]==highest_degree_vert){
    //             siblings[verts_to_erase].erase(siblings[verts_to_erase].begin() +i);
    //         }
    //     }
    // }

    //remove u and v info from all other vertex
    for (int verts_to_erase: siblings[v]) {
      for (int i = 0; i < siblings[verts_to_erase].size(); i++) { //go thru the verts that need to forget that v exist
        if (siblings[verts_to_erase][i] == v) {
          siblings[verts_to_erase].erase(siblings[verts_to_erase].begin() + i);
        }
      }
    }
    for (int verts_to_erase: siblings[u]) {
      for (int i = 0; i < siblings[verts_to_erase].size(); i++) { //go thru the verts that need to forget that u exist
        if (siblings[verts_to_erase][i] == u) {
          siblings[verts_to_erase].erase(siblings[verts_to_erase].begin() + i);
        }
      }
    }
    siblings[u].clear(); //clear all neighbor info from vertex u
    siblings[v].clear(); //clear all neighbor info from vertex u

  }

  outVC2 << "APPROX-VC2: ";
  // for (int k : assignments){
  //     std::cout<<k<<",";
  // }

  std::sort(assignments.begin(), assignments.end()); //sort the sat assignments from lowest to highest

  for (int i = 0; i < assignments.size(); i++) { //print out the solution
    if (i == assignments.size() - 1) {
      outVC2 << assignments[i] << std::endl;
    } else {
      outVC2 << assignments[i] << ",";
    }

  }
  // std::cout<<outVC2.str()<<std::flush;
  return assignments.size();
}

//approx vc2 thread time
void * APPROX_VC_2Thread(void * arg) {
  clock_t ON, OFF;
  double duration;
  // std::cout<<pthread_self()<<std::endl<<std::flush;
  int test_v = * (int * ) arg;
  // while (!g_cancel){
  //     if (g_RunOnce){
  ON = clock();
  // std::cout<< "inside vc2thread"<< v<<std::endl<<std::flush;
  APPROX_VC_2(std::ref(undirected_graph -> adj_l), test_v);
  OFF = clock();

  duration = (double)(OFF - ON) / CLOCKS_PER_SEC;
  outVC2_time << "vc2 time: " << duration;
  // g_cancel=true;
  //     }
  // }
  // std::cout<<"vc2 thread exit..."<<std::endl;
  //cout << "approxVCThread 2----------" << duration << endl;
  return nullptr;

}

void * threadInOut(void * somearg) {

  std::stringstream out;
  std::ofstream myfile;

//   int Optimal_CoverSize, Approx1_CoverSize, Approx2_CoverSize;
  //read from stdin until EOF
  while (!std::cin.eof()) {
    std::string line;
    std::getline(std::cin, line);

    if (line.empty()) {
      // std::cout<<"donnneee"<<std::endl;
      g_cancel = true;
      continue;
    }
    try {
      int edge_count = 0;

      if (line.find('V') == 0) {
        v = std::stoi(line.substr(2));
        delete undirected_graph;
        undirected_graph = new create_graph(v);
      }
      if (line.find('E') == 0) {
        std::size_t spec_char = line.find('<');
        int flag = 0;
        while (spec_char != std::string::npos) {
          std::size_t spec_char1 = line.find('>', spec_char);
          std::string edge_input = line.substr(spec_char + 1, spec_char1 - spec_char - 1);
          std::size_t seperator = edge_input.find(',');
          int edge_0 = std::stoi(edge_input.substr(0, seperator));
          int edge_1 = std::stoi(edge_input.substr(seperator + 1));
          if (edge_0 > v or edge_1 > v or edge_0 < 0 or edge_1 < 0 or edge_0 == edge_1) {
            std::cout << "Error: Incorrect edge input." << std::endl;
            flag = 1;
            break;
          }
          spec_char = line.find('<', spec_char + 1);
        }
        edge_count = 0;
        if (flag == 0) {
          std::size_t spec_char2 = line.find('<');
          while (spec_char2 != std::string::npos) {
            std::size_t spec_char3 = line.find('>', spec_char2);
            std::string edge_input = line.substr(spec_char2 + 1, spec_char3 - spec_char2 - 1);
            std::size_t seperator = edge_input.find(',');
            int edge_0 = std::stoi(edge_input.substr(0, seperator));
            int edge_1 = std::stoi(edge_input.substr(seperator + 1));
            edge_add(undirected_graph -> adj_l, edge_0, edge_1);
            edge_count++;
            spec_char2 = line.find('<', spec_char2 + 1);
          }
        }
        g_cancel = 0;
        g_RunOnce = 1;
        //thread creation 
        int c;

        pthread_t thread1, thread2, thread3, thread4;

        int test_v = v;
        // std::cout<< "inside inouthtread"<< v<<std::flush;
        c = pthread_create( & thread1, NULL, & solveVertexCoverThread, & test_v);
        c = pthread_create( & thread2, NULL, & APPROX_VC_1Thread, & test_v);
        c = pthread_create( & thread3, NULL, & APPROX_VC_2Thread, & test_v);

        pthread_create( & thread4, NULL, & watchdog, & thread1);

        pthread_join(thread1, NULL);
        // std::cout<<"colleted cnf thread..."<<std::endl;
        pthread_join(thread2, NULL);
        // std::cout<<"colleted approx1 thread..."<<std::endl;
        pthread_join(thread3, NULL);
        // std::cout<<"colleted approx2 thread..."<<std::endl;

        g_cancel = 1;
        pthread_join(thread4, NULL);
        // std::cout<<"colleted watchdog thread..."<<std::endl;

        std::cout << outVC.str() << std::flush; //read the string streams created in the threads
        std::cout << outVC1.str() << std::flush;
        std::cout << outVC2.str() << std::flush;

        outVC.flush();
        outVC1.flush();
        outVC2.flush();

        // Optimal_CoverSize= solveVertexCover(v, undirected_graph -> adj_l, edge_count);

        // //approx_vc1 here?
        // Approx1_CoverSize=APPROX_VC_1(undirected_graph -> adj_l,v );

        // out<<"Appox VC1 ratio: "<<(float)Approx1_CoverSize/(float)Optimal_CoverSize<<std::endl;

        // Approx2_CoverSize=APPROX_VC_2(undirected_graph -> adj_l,v );
        // out<<"Appox VC2 ratio: "<<(float) Approx2_CoverSize/(float) Optimal_CoverSize<<std::endl;
        // // std::cout<<std::cin.eof()<<"----217---"<<std::endl;

        //============SAVE RESULTS TO FILES================= //

        // myfile.open("size_"+ std::to_string(v) + "output.txt", std::ios_base::app); //append to the end of this file
        // myfile<<out.rdbuf(); //read the buffer
        // out.flush(); 
        // myfile.close();
        myfile.open("size_" + std::to_string(v) + "_TIMEoutput.txt", std::ios_base::app); //append to the end of this file
        myfile << outVC_time.str() << " " << outVC1_time.str() << " " << outVC2_time.str() << std::endl << std::flush; //read the buffer
        myfile.close();
        outVC_time.flush();
        outVC1_time.flush();
        outVC2_time.flush();
        // myfile.open("size_"+ std::to_string(v) + "_TIMEoutput.txt", std::ios_base::app); //append to the end of this file
        // myfile<<out.rdbuf(); //read the buffer
        // out.flush(); 
        // myfile.close();
        // myfile.open("size_"+ std::to_string(v) + "_TIMEoutput.txt", std::ios_base::app); //append to the end of this file
        // myfile<<out.rdbuf(); //read the buffer
        // out.flush(); 
        // myfile.close();

      }

    } catch (const std::exception & e) {
      std::cout << "Error: Invalid Input" << '\n';
    }
  }
  return nullptr;
}

int main(int argc, char ** argv) {

  int c;
  pthread_t threadinout;
  // g_cancel=false;

  c = pthread_create( & threadinout, NULL, threadInOut, NULL);
  pthread_join(threadinout, NULL);
  return 0;
}