#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string.h>
#include <queue>
#include <algorithm>
#include <memory>
#include "minisat/core/SolverTypes.h"
#include "minisat/core/Solver.h"
#include <pthread.h>
#include <time.h>

using namespace std;
using namespace Minisat;

#define MAX 100000

#define INFINITY 200000

vector<int> find_vertexcover(int, vector<int>, int);
void CNF_SAT_VC(int, vector<int>);
void APPROCX_VC_1(int, vector<int>);
void APPROCX_VC_2(int, vector<int>);

//int pthread_getcpuclockid(pthread_t thread,clockid_t *clock_id);
vector<int> vec_V;
vector<int> vec_E;
int vertex_num;

vector<int> output_0;
vector<int> output_1;
vector<int> output_2;


bool is_element_in_vector(vector<int> v, int element) {
	vector<int>::iterator it;

	it = find(v.begin(), v.end(), element);
	if (it != v.end()) {
		return true;
	}

	else {
		return false;
	}
}

void* threadAPPROX_VC_1(void* arg) {
	//cout << "1" << endl;

	clock_t start, finish;
	double duration;

	start = clock();
	APPROCX_VC_1(std::ref(vertex_num), std::ref(vec_E));
	finish = clock();

	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "the duration time is: " << duration << endl;

}


void APPROCX_VC_1(int vertex_num, vector<int>edge) {
	int cnt[vertex_num];
	int edge_num = edge.size();
	int max = 0;
	std::vector<int> new_edge = edge;
	//std::vector<int> output_1;

	while (!edge.empty()) {
		for (int i = 0; i < vertex_num; i++) {
			cnt[i] = 0;
		}
		for (int i = 0; i < edge.size(); i++) {
			if (edge[i] != INFINITY)
				cnt[edge[i]]++;

		}

		max = max_element(cnt, cnt + vertex_num) - cnt;
		output_1.push_back(max);
		new_edge.clear();


		for (int i = 0; i < edge.size(); i++) {
			if (edge[i] == max || edge[i + 1] == max) {
				//continue;
				edge[i] = INFINITY;
				edge[i + 1] = INFINITY;
				//edge_num -= 2;
			}
			else {
				new_edge.push_back(edge[i]);
				new_edge.push_back(edge[i + 1]);
			}
			i++;
		}
		edge = new_edge;
	}

	sort(output_1.begin(), output_1.end());
}

void* threadAPPROX_VC_2(void* arg) {
	clock_t start, finish;
	double duration;

	start = clock();
	APPROCX_VC_2(ref(vertex_num), ref(vec_E));
	finish = clock();

	duration = (double)(finish - start)/CLOCKS_PER_SEC;
	//cout << "the duration time is: " << duration << endl;

}

void APPROCX_VC_2(int vertex_num, vector<int>edge) {

	//std::vector<int> output_2;
	vector<int> new_edge = edge;


	while (!edge.empty()) {

		output_2.push_back(edge[0]);
		output_2.push_back(edge[1]);

		new_edge.clear();
		for (unsigned i = 0; i < edge.size() - 1; ++i) {
			if (i % 2 == 0) {
				if (edge[i] == edge[0] || edge[i + 1] == edge[1] || edge[i] == edge[1] || edge[i + 1] == edge[0]) {
					continue;
				}
				else {
					new_edge.push_back(edge[i]);
					new_edge.push_back(edge[i + 1]);
				}
			}
		}
		edge = new_edge;
	}
	sort(output_2.begin(), output_2.end());
}

void* threadCNF_SAT_VC(void* arg) {

	//CNF_SAT_VC(std::ref(vertex_num), std::ref(vec_E));
	clock_t start, finish;
	double duration;

	start = clock();
	CNF_SAT_VC(std::ref(vertex_num), std::ref(vec_E));
	finish = clock();

	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "the duration time is: " << duration << endl;

}


void CNF_SAT_VC(int vertex_num, vector<int>edge) {
	Solver solver;

	vector<int> result;
	//vector<int> output_0;

	int num = vertex_num;

	vector<int> vec_edge = edge;

	int hi = vertex_num;
	int low = 1;
	//int mid;
	vector<int> tmp = { -1 };

	while (hi >= low) {
		int mid = (hi + low) / 2;

		result = find_vertexcover(num, vec_edge, mid);

		bool nosat = std::equal(result.begin(), result.end(), tmp.begin());
		if (not nosat) {
			hi = mid - 1;
			output_0.clear();
			output_0 = result;
		}
		else {
			low = mid + 1;
		}
	}
	if (vec_edge.empty()) {
		output_0.clear();
	}
	sort(output_0.begin(), output_0.end());

}
std::vector<int> find_vertexcover(int vertex_num, vector<int>edge, int k) {

	Solver solver;

	Var propositions[vertex_num][k];

	//vector<vector<Minisat::Lit> > vec;

	for (int i = 0; i < vertex_num; i++) {
		for (int j = 0; j < k; j++) {
			propositions[i][j] = solver.newVar();
		}
	}

	// clause1: at least (exactly only one) one vertex is the ith vertex in the vertex cover, i in [1,k]

	for (int i = 0; i < k; i++) {
		vec<Lit> clause1;

		for (int j = 0; j < vertex_num; j++) {
			clause1.push(mkLit(propositions[j][i]));
		}

		solver.addClause(clause1);
		clause1.clear();

	}

	// clause 2. No one vertex can appear twice in a vertex cover

	for (int m = 0; m < vertex_num; m++) {
		for (int p = 0; p < k - 1; p++) {
			for (int q = p + 1; q < k; q++) {
				solver.addClause(mkLit(propositions[m][p], true), mkLit(propositions[m][q], true));
			}
		}
	}
	// clause 3. No more than one vertex appears in the mth position of the vertex cover.

	for (int m = 0; m < k; m++) {
		for (int p = 0; p < vertex_num - 1; p++) {
			for (int q = p + 1; q < vertex_num; q++) {
				solver.addClause(mkLit(propositions[p][m], true), mkLit(propositions[q][m], true));
			}
		}
	}

	// 4. Every edge is incident to at least one vertex in the vertex cover
	for (int i = 0; i < edge.size(); i++) {
		vec<Lit> clause4;
		for (int m = 0; m < k; m++) {
			clause4.push(mkLit(propositions[edge[i]][m]));

			clause4.push(mkLit(propositions[edge[i + 1]][m]));
		}
		solver.addClause(clause4);
		clause4.clear();
		i++;
	}

	std::vector<int> vertex_cover;


	if (solver.solve()) {

		for (int i = 0; i < vertex_num; i++) {
			for (int j = 0; j < k; j++) {
				if (toInt(solver.modelValue(propositions[i][j])) == 0) {
					vertex_cover.push_back(i);
				}
			}
		}

		return vertex_cover;
		//sort	
	}
	else {
		return { -1 };
	}
}

void output() {
	double ratio_0, ratio_1, ratio_2;

	cout << "CNF-SAT-VC: ";
	for (int i = 0; i < output_0.size(); ++i) {
		if (i < output_0.size()-1) {
			cout << output_0[i] << ",";
		}
		else {
			cout << output_0[i];
		}
	}
	cout << "\n";

	cout << "APPROX-VC-1: ";
	for (int i = 0; i < output_1.size(); ++i) {
		if (i < output_1.size()-1) {
			cout << output_1[i] << ",";
		}
		else {
			cout << output_1[i];
		}		
	}
	cout << "\n";

	cout << "APPROX-VC-2: ";
	for (int i = 0; i < output_2.size(); ++i) {
		if (i < output_2.size()-1) {
			cout << output_2[i] << ",";
		}
		else {
			cout << output_2[i];
		}
	}
	cout << "\n";

	output_0.clear();
	output_1.clear();
	output_2.clear();
}

void* threadIO(void* arg) {
	//cout << "io begin" << endl;

	while (!::cin.eof()) {

		string line;

		while (getline(cin, line)) {
			char arr[MAX];

			strcpy(arr, line.c_str());

			if (arr[0] == 'V') {
				//cout << "vertex\n";
				//cout << line << endl; //a3: edit it to print the graph that it has read before accepting an s command.

				int num = 0;
				for (int i = 0; arr[i] != '\0'; ++i) {

					if (arr[i] >= '0' && arr[i] <= '9')

						num = num * 10 + arr[i] - '0';
				}


				for (int j = 0; j < num; ++j) {
					vec_V.push_back(j);
				}
				//for (vector<int>::iterator it = vec_V.begin(); it != vec_V.end(); it++) {
				//	//cout << *it << " ";
				//}//print vector
				while (vec_V.empty()) {
					//cout << "Error: NO vertex is enterd." << endl;
					break;
				}
				break;
			}


			else if (arr[0] == 'E') {
				int num = 0;
				int d_edge = 0;
				//cout << line << endl;//a3: edit it to print the graph that it has read before accepting an s command.
				for (int i = 0; arr[i] != '\0'; i++) {

					if (arr[i] > '0' && arr[i] <= '9') {
						/*cout << arr[i] << "\n";*/
						num = num * 10 + arr[i] - '0';
					}
					else {
						if (arr[i] == '0' && num == 0) {
							vec_E.push_back(0);
						}
						if (arr[i] == '0' && (num != 0)) {
							num = num * 10 + arr[i] - '0';
							vec_E.push_back(num);
							num = 0;
						}
						if (arr[i] != '0' && (num != 0)) {
							vec_E.push_back(num);
							num = 0;
						}

					}
				}
				// for (int i = 0; i < vec_E.size(); i++) {

				// 	if (is_element_in_vector(vec_V, vec_E[i]) == false) {
				// 		//cout << vec_E[i];
				// 		cout << "Error:some vertex in  edge not exists.";// check the error of edge: not a vertex
				// 		vec_E.erase(vec_E.begin() + i);
				// 		return 0;
				// 	}
				// }
				while (vec_E.empty()) {
					//cout << "Error: NO edge is enterd." << endl;
					//return 0;
					break;
				}
				vertex_num = vec_V.size();
				
				int s;
				pthread_t threadio, thread0, thread1, thread2;

				s = pthread_create(&thread0, NULL, threadCNF_SAT_VC, NULL);
				s = pthread_create(&thread1, NULL, threadAPPROX_VC_1, NULL);
				s = pthread_create(&thread2, NULL, threadAPPROX_VC_2, NULL);

				pthread_join(thread0, NULL);
				pthread_join(thread1, NULL);
				pthread_join(thread2, NULL);

				output();

				vec_E.clear();
				vec_V.clear();

			}
		}
	}
}


int main(int argc, char** argv) {

	int s;
	pthread_t threadio, thread0, thread1, thread2;

	//clockid_t cid0,cid1,cid2;

	//vector<int> result_cnf, result_a1, result_a2;

	s = pthread_create(&threadio, NULL, threadIO, NULL);
	//s = pthread_create(&thread0, NULL, threadCNF_SAT_VC, NULL);
	//s = pthread_create(&thread1, NULL, threadAPPROX_VC_1, NULL);
	//s = pthread_create(&thread2, NULL, threadAPPROX_VC_2, NULL);


	pthread_join(threadio, NULL);
	//pthread_join(thread0, NULL);
	//pthread_join(thread1, NULL);
	//pthread_join(thread2, NULL);
	return 0;
}
