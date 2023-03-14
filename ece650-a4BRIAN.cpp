// Compile with c++ ece650-a2cpp -std=c++11 -o ece650-a2
#include <iostream>
#include <algorithm>
#include <sstream>
#include <vector>
#include <queue>
#include <stack>

// defined std::unique_ptr
#include <memory>
// defines Var and Lit
#include "minisat/core/SolverTypes.h"
// defines Solver
#include "minisat/core/Solver.h"

std::vector<std::string> seperate_str(std::string s, std::string delimiter=" ") {
    std::vector<std::string> v={};

	size_t pos = 0;
	size_t prev_pos = 0;
	std::string token;
	while ((pos = s.find(delimiter,prev_pos)) != std::string::npos) {
    	token = s.substr(prev_pos, pos-prev_pos);
    	prev_pos=pos+delimiter.length(); //update position for beginning of next section to check
        v.push_back(token); //save my token to end of vector
    }
	token=s.substr(prev_pos);
    v.push_back(token);
    return v;
    
}

void strip(std::string &s, std::string delimiter=" "){
    size_t pos = 0;
	size_t prev_pos = 0;
    while ((pos = s.find(delimiter,prev_pos)) != std::string::npos) {
        s.erase(pos,1);
        prev_pos=pos+delimiter.length();
    }
}
void PrintSiblingDatabase(std::vector<std::vector<int>> siblings){
    for( int i=0; i<siblings.size();i++){
        std::cout<<"Vert number: "<<i<<std::endl;
        for (int n: siblings[i]){
            std::cout<<n<<std::endl;
        }
        
    }
}
void ClearSiblingDatabase(std::vector<std::vector<int>> &siblings){
    for( int i=0; i<siblings.size();i++){
        // std::cout<<"Vert number: "<<i<<std::endl;
        siblings[i].clear();
        
    }
    siblings.clear();
    // std::cout<<"After clear fn: "<<siblings.empty()<<std::endl;
}

void parse_command(std::string line, char &cmd, std::vector<int> &coords, int &start_vert, int &end_vert, int &number_of_verts,std::vector<std::vector<int>> siblings, std::string &errmsg){ //breakdown the given user command and update info in variables
    std::vector<std::string> command,coords_str;
    std::string coords_RAW;

    command=seperate_str(line); //return vector holding tokens created with " " delim
    cmd=command[0][0]; //save the command type
        if (command[0]=="E") {
            // std::cout<<"E command"<<std::endl;
            coords.clear();
            if (number_of_verts>=1){
                if (command[1].length()>3){
                    // std::cout<<command[1].length()<<std::endl;
                    coords_RAW=command[1].erase(0,2); //REMOVE "{<..."
                    coords_RAW.erase(coords_RAW.length()-2,2); //REMOVE "...>}"
                    // std::cout<< coords_RAW<<std::endl;
                    strip(coords_RAW, "<"); //REMOVE ALL INSTANCE OF < AND >
                    strip(coords_RAW, ">");
                    // std::cout<< coords_RAW<<std::endl;
                    coords_str=seperate_str(coords_RAW,","); //coords_str holds the vertex number as VECTORstr

                    for(std::string n: coords_str){
                        coords.push_back(std::stoi(n)); //coords holds the vertex number as VECTORint
                        // std::cout<<"coords holds:"<<std::stoi(n)<<"\n";
                    }

                    //check if coords allowed
                    for(int n: coords){
                        //    std::cout<<n<<std::endl; //print contents of the coords vector
                        if (!(n<=number_of_verts && n>0 )){
                            errmsg="vertex referenced not allowed";
                            // std::cout<<n<<number_of_verts<<std::endl;
                            ClearSiblingDatabase(siblings);
                            // std::cout<<"siblings cleared"<<std::endl;
                            // std::cout<<"After clear fn: "<<siblings.empty()<<std::endl;
                             //prevent s cmd call by invalidating prev info
                            throw 99;
                        }
                    }    
                } else {
                    // std::cout<<command[1]<<std::endl;
                    coords.clear();
                    coords.resize(1);
                }
            } else {
                errmsg="V not called yet";
                throw 99;
            }



        } else if (command[0]=="s" ){
            // std::cout<<"S command"<<std::endl;
            start_vert=std::stoi(command[1]);
            end_vert=std::stoi(command[2]);
            // for(int n: coords){
            //    std::cout<<n<<std::endl; //print contents of the coords vector
            //    if (!(n<=number_of_verts)){
            //        errmsg="vertex referenced not allowed";
            //         throw 99;
            //    }
            // }
            // std::cout<<"in parse: "<<siblings.empty()<<std::endl;
            // PrintSiblingDatabase(siblings);
            if (!siblings.empty() && (siblings.size()!=0) ){
                if ( (siblings[start_vert].empty() || siblings[end_vert].empty() ) && (start_vert!=end_vert)){ //case of s 1 2
                    errmsg="start or end vert not in tree";
                    throw 99;
                } else if ( (start_vert!=end_vert) ){
                    if ( !(start_vert<= number_of_verts && start_vert>0) ){
                        errmsg="start vert not allowed";
                        throw 99; 
                    }
                    if ( !(end_vert<= number_of_verts && end_vert>0) ){
                        errmsg="end vert not allowed";
                        throw 99; 
                    }
                } else if ( start_vert==end_vert ){ //case of s 1 1 check if allowed
                    //check if coords allowed
                    if ( !(start_vert<= number_of_verts && start_vert>0) ){
                        errmsg="start and end vert not allowed";
                        throw 99; 
                    }
                }
            } else {
                errmsg="no edges specified";
                throw 99;
            }

            // std::cout<<start_vert<<":"<<end_vert<<std::endl;
        }else if (command[0]=="V"){
            // std::cout<<"V command"<<std::endl;
            if ( std::stoi(command[1])>1 ){
                number_of_verts= std::stoi(command[1]);
                siblings.clear();
            } else {
                errmsg="Vert must be >1";
                siblings.clear();
                number_of_verts=0;
                coords.clear();
                coords.resize(1);
                //prevent s cmd call by invalidating prev info

                throw 99;
            }
            // std::cout<<number_of_verts<<std::endl;
        }else{
            errmsg="Command>"+command[0]+"<not recognized";
            throw 99;
        }

}

void SiblingDatabase(std::vector<int> coords,std::vector<std::vector<int>> &siblings){
    // std::cout<< "hi"<<std::endl;
    
    if (coords.size()>1) {
        for (int i=0; i<coords.size(); i+=2){
            //assuming not repeat edges or bc unidirected graph, no need to check if already exist in vector
            siblings.at(coords[i]).push_back(coords[i+1]); //{{ },...{i+1},...{}...} value at position
                                                        //{{0},...{i},.....{}...} position in vector
            siblings.at(coords[i+1]).push_back(coords[i]); //ADD THE REVERSE OF THE COORD AS SIBLINGS
            // std::cout<<coords[i]<<":"<<coords[i+1]<<std::endl;
        }
    }
}

void APPROX_VC_1(std::vector<std::vector<int>> siblings){
    // Taking the highest degree vertex and removing children
    int highest_degree=0;
    int highest_degree_vert=0;
    std::vector<int> assignments;
    bool sibling_info_empty=false;

    
    while(!sibling_info_empty){
        // PrintSiblingDatabase(siblings);
        for(std::vector<int> vert: siblings){
            if(vert.size()==0){
                sibling_info_empty=true;
                
            }else{
                sibling_info_empty=false;
                break;
            }
        }

        if (sibling_info_empty==true){
            break;
        }

        for (int vert=1; vert<siblings.size(); vert++){
            if (siblings[vert].size()>0){
                if (siblings[vert].size() >highest_degree){
                    highest_degree=siblings[vert].size();
                    highest_degree_vert=vert; //update the highest degree vert 
                } 
            }
        }

        assignments.push_back(highest_degree_vert);
        for(int verts_to_erase :siblings[highest_degree_vert]){ //remove highest_degree vertex from all othe rvertex's info
            for(int i=0; i<siblings[verts_to_erase].size(); i++){
                if (siblings[verts_to_erase][i]==highest_degree_vert){
                    siblings[verts_to_erase].erase(siblings[verts_to_erase].begin() +i);
                }
            }
        }
        siblings[highest_degree_vert].clear();
        highest_degree=0;
        
    }

    
    
        std::cout<<"PRINTING VC1..."<<std::endl;
        // for (int k : assignments){
        //     std::cout<<k<<",";
        // }

        std::sort(assignments.begin(),assignments.end()); //sort the sat assignments from lowest to highest

        for(int i=0; i<assignments.size(); i++){ //print out the solution
            if (i==assignments.size()-1){
                std::cout<<assignments[i]<<std::endl;  
            } else {
                std::cout<<assignments[i]<<" ";  
            }
            
        }

}

void APPROX_VC_2(std::vector<std::vector<int>> siblings){
    // Taking the highest degree vertex and removing children
    int highest_degree=0;
    int highest_degree_vert=0;
    std::vector<int> assignments;
    bool sibling_info_empty=false;

    int u,v;


    
    while(!sibling_info_empty){
        // PrintSiblingDatabase(siblings);
        for(std::vector<int> vert: siblings){
            if(vert.size()==0){
                sibling_info_empty=true; //checking if all siblings are covered
                
            }else{
                sibling_info_empty=false;
                break;
            }
        }

        if (sibling_info_empty==true){
            break;
        }

        //pick an edge to use

        for (int vert=1; vert<siblings.size(); vert++){
            if (siblings[vert].size()>0){
                u=vert;

                v=siblings[vert].back();
                siblings[vert].pop_back();
                break;
            }
        }

        assignments.push_back(v);
        assignments.push_back(u);
        std::cout<<u<<" : "<<v<<std::endl;
        
        // for(int verts_to_erase :siblings[highest_degree_vert]){ //remove highest_degree vertex from all othe rvertex's info
        //     for(int i=0; i<siblings[verts_to_erase].size(); i++){
        //         if (siblings[verts_to_erase][i]==highest_degree_vert){
        //             siblings[verts_to_erase].erase(siblings[verts_to_erase].begin() +i);
        //         }
        //     }
        // }

        //remove u and v info from all other vertex
        for(int verts_to_erase :siblings[v]){ 
            for(int i=0; i<siblings[verts_to_erase].size(); i++){ //go thru the verts that need to forget that v exist
                if (siblings[verts_to_erase][i]==v){
                    siblings[verts_to_erase].erase(siblings[verts_to_erase].begin() +i);
                }
            }
        }
        for(int verts_to_erase :siblings[u]){ 
            for(int i=0; i<siblings[verts_to_erase].size(); i++){ //go thru the verts that need to forget that u exist
                if (siblings[verts_to_erase][i]==u){
                    siblings[verts_to_erase].erase(siblings[verts_to_erase].begin() +i);
                }
            }
        }
        siblings[u].clear();//clear all neighbor info from vertex u
        siblings[v].clear();//clear all neighbor info from vertex u
        
    }

    
    
        std::cout<<"PRINTING VC2..."<<std::endl;
        // for (int k : assignments){
        //     std::cout<<k<<",";
        // }

        std::sort(assignments.begin(),assignments.end()); //sort the sat assignments from lowest to highest

        for(int i=0; i<assignments.size(); i++){ //print out the solution
            if (i==assignments.size()-1){
                std::cout<<assignments[i]<<std::endl;  
            } else {
                std::cout<<assignments[i]<<" ";  
            }
            
        }

}

std::vector<int> solve(int start_vert, int end_vert ,std::vector<std::vector<int>> siblings,int number_of_verts){

    std::queue<int> MyQueue;
    std::vector<int> prev;
    std::vector<bool> visted;
    visted.resize(number_of_verts+1,false);//index 0not used. set entire vector to false
    prev.resize(number_of_verts+1);
    int node;

    visted[start_vert]=true;
    MyQueue.push(start_vert);
    while ( !MyQueue.empty() ) {
        node= MyQueue.front();
         MyQueue.pop();
        //  std::cout<<node<<std::endl;
        // std::cout<<"START"<<std::endl;
        // for(int n: siblings[node]){
        //        std::cout<<n<<std::endl;
        //    }
        //    std::cout<<"END"<<std::endl;


        for(int next:siblings[node]){
            
            if (!visted[next]){
                // std::cout<<next<<std::endl;
                MyQueue.push(next);
                visted[next]=true;
                prev[next]=node;
            }
        }
        // std::cout<<"BRIANN"<<std::endl;
    }
    return prev;
}
void reverseQueue(std::queue<int> &Queue){
    // https://www.geeksforgeeks.org/reversing-a-queue/
    std::stack<int> Stack;
    //3-point turn with queue as a car :)
    while (!Queue.empty()) {
        Stack.push(Queue.front());
        Queue.pop();
    }
    while (!Stack.empty()) {
        Queue.push(Stack.top());
        Stack.pop();
    }
}

std::queue<int> ReconstructPath(int start_vert, int end_vert, std::vector<int> prev, std::string &errmsg){
    std::queue<int> path;

    int n= end_vert;
    while( n!=0){
        path.push(n);
        n=prev[n];
    }

    // while(!path.empty()){
    //     std::cout<< path.front()<<std::endl;
    //     path.pop();
    // }

    reverseQueue(path);
    // while(!path.empty()){
    //     std::cout<< path.front()<<std::endl;
    //     path.pop();
    // }

    if (path.front()==start_vert){ //path is in reverse order still ie.) [end,...,start]
        return path;
    } else{
        while(!path.empty()){
        path.pop(); //empty queue
        }
        errmsg="No path exist";
        throw 99;
        return path; //return empty path
    }
}

int main(int argc, char** argv) {

    // read from stdin until EOF
    int start_vert=0, end_vert=0, number_of_verts=0;
    std::vector<std::vector<int>> siblings;
    std::vector<int> prev;
    std::queue<int> path;
    std::string errmsg =" ";
    std::vector<int> coords;
    char cmd;

    while (!std::cin.eof()) {
        try {
            // read a line of input until EOL and store in a string
            std::string line;
            std::getline(std::cin, line);
            if (line.length()==0){ //exit program if empty input
                return 1;
            }
            // create an input stream based on the line 
            // we will use the input stream to parse the line
            std::istringstream input(line);

            // we expect each line to contain a list of numbers
            // this vector will store the numbers.
            // they are assumed to be unsigned (i.e., positive)
            
            // std::cout<<"read:"<<line<<std::endl;
            // std::cout<<start_vert<<std::endl;
            // std::cout<<end_vert<<std::endl;
            // std::cout<<number_of_verts<<std::endl;
            // std::cout<<cmd<<std::endl;
            parse_command(line,cmd,coords,start_vert,end_vert,number_of_verts, siblings, errmsg);
            //plus 1 bc index 0 wont be used

            // siblings.resize(number_of_verts+1);// PLUS 1 BC INDEX 0 NEVER USED

            // std::cout<<"SIZE"<<siblings.size()<<std::endl;
            // for( int i=0; i<siblings.size();i++){
            //             for (int n: siblings[i]){
            //                 std::cout<<n<<std::endl;
            //             }
            //             std::cout<<"next"<<std::endl;
            //         }

            // for( int i=0; i<siblings.size();i++){
            //     siblings[i].resize(10);
            // }
            // std::cout<<"{"<<std::endl;
            // std::cout<<start_vert<<std::endl;
            // std::cout<<end_vert<<std::endl;
            // std::cout<<number_of_verts<<std::endl;
            // std::cout<<cmd<<std::endl;
            // std::cout<<"}"<<std::endl;
            // std::cout<<"in main,after parse smd: "<<siblings.empty()<<std::endl;
            
            if (cmd== 'V'){
                
                siblings.resize(number_of_verts+1);
                siblings.clear();

                // std::cout<<"V CMD1"<<std::endl;
            } else if (cmd== 'E'){
            //     for(int n: coords){
            //        std::cout<<"COORD: "<<n<<std::endl;
            //    }
                if (coords.size()==1){
                    std::cout<<"NULL COVER"<<std::endl;
                } else {
                    siblings.resize(number_of_verts+1);
                    SiblingDatabase(coords,siblings);
                    // std::cout<< "bye"<<std::endl;
                    // PrintSiblingDatabase(siblings);
                    APPROX_VC_1(siblings);
                    APPROX_VC_2(siblings);
                }
                
            } else if (cmd=='s'){
                // std::cout<<start_vert<<std::endl;
                // std::cout<<end_vert<<std::endl;
                // std::cout<<number_of_verts<<std::endl;
                // std::cout<<cmd<<std::endl;
                // PrintSiblingDatabase(siblings);
                // std::cout<< "bye"<<std::endl;
                prev= solve(start_vert,end_vert,siblings,number_of_verts);
                // std::cout<<"START"<<std::endl;
                // for(int n: prev){
                //     std::cout<<n<<std::endl;
                // }
                // std::cout<<"END"<<std::endl;
                path=ReconstructPath(start_vert,end_vert,prev,errmsg);
                while(!path.empty()){
                    if (path.size()!=1){
                        std::cout<< path.front()<<"-";
                    } else {
                        std::cout<< path.front()<<std::endl;
                    }
                    
                    path.pop();
                }
            }
            // std::cout<<"FINAL line: "<<siblings.empty()<<std::endl;
        
        } 
        catch(const std::exception& ex){
                std::cout<<"Error: "<< ex.what()<<"\n";
            }
        catch(...){
                std::cout<<"Error: CustomMsg->"<< errmsg<<"\n";
            }
        
    }
    
}
