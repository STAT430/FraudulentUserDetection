
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <list>
#include <vector>
#include <assert.h>
#include <string.h>
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <pthread.h>
#include <cmath>
#include <time.h>

#define GCC_VERSION (__GNUC__ * 10000 \
+ __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

#if GCC_VERSION >= 40300
#include <tr1/unordered_map>
using namespace std::tr1;
//namespace std {
//    std::string to_string(size_t n) {
//        std::ostringstream s;
//        s << n;
//        return s.str();
//    }
//}
#define hash_map unordered_map

#else
#include <unordered_map>
#endif
using namespace std;

// random generator function:
ptrdiff_t myrandom (ptrdiff_t i) { return rand()%i;}

// pointer object to it:
ptrdiff_t (*p_myrandom)(ptrdiff_t) = myrandom;

long randint(long min, long max){

    return long(rand() / (RAND_MAX + 0.0) * (max - min) + min);

}

bool my_comparator ( const pair<double, int>& l, const pair<double, int>& r)
{
    return l.first > r.first;
}

class Data{
public:
    typedef unsigned long int vertex;

    long N;

    //adjacency list
    //weighted graph
    unordered_map<vertex,set<pair<vertex, double> > > network_map;

    //1 : weighted, 0: unweighted
    int weighted_graph;

    //priors
    //probability of being benign for each node
    double *prior;

    //network file
    char *network_file;

    char *post_file;
    char *prior_file;
    
    //char *post_file_iter;

    //during computation, post is treated as the normalized multiplication of incoming messages of each node.
    //
    double *post;
    double *post_pre;

    //to support different input format
    char *train_set_file;
    
    double alpha;

    //parameters of PMR model when input is a set of labeled nodes, instead of priors
    double theta_pos;
    double theta_neg;
    double theta_unl;

    //this weight is used when the input is an unweighted graph
    double weight;

    double max_iter;

    int* ordering_array;

    int num_threads;

	set<vertex> neg_train_set;
    // directed graph edges
    //int** graph_array;
    //unordered_map<pair<vertex, vertex>, pair<int, double>, pairhash > graph_array;
    
    set<long> valid_user;
    //edge types
    //1 : ourgoing
    //2 : incoming
    //3 : bidirectional

    class CIA_arg{
    public:
        Data * data_pointer;
        int current_thread;
        CIA_arg(){}
    };

    Data(){


    }

    void add_edge(vertex node1, vertex node2, double w){
        // add edge (node1, node2)
        
        // no self loops
        if(node1==node2){
            return;
        }
		//add node2 to the adjacency list of node1
		network_map[node1].insert(make_pair(node2, w));
        
    }
    
    /* Read in the social graph */
    //the format for the social graph is
    //each line corresponds to an edge, e.g, 3 2 0.8
    //each edge in the graph appears twice, e.g.,
    //3 2 0.8
    //2 3 0.9
    void read_network(){
        
        ifstream in(network_file,ifstream::in);
        assert(in);
        
        string line;
        vertex node1,node2, max_node=0;
        double w;
        
        //read edges
        while(getline(in,line)!=NULL){
            node1=(vertex)atol(strtok((char *)line.c_str()," \n\t\r"));
            node2=(vertex)atol(strtok(NULL," \n\t\r"));
            if (weighted_graph == 1) {
                w=(double)atof(strtok(NULL," \n\t\r"));
            }
            else{
                w = weight;
            }
            //cout << node1 << " " << node2 << " " << w << endl;
            add_edge(node1, node2,w);
            if(node1 > max_node){
                max_node = node1;
            }
            if(node2 > max_node){
                max_node = node2;
            }
            
        }
        
        in.close();
        //number of nodes in the graph
        N = network_map.size();
        //N = max_node + 1;
        //cout<<N<<endl;
        
        //allocate space for final scores
        post = (double*) malloc(sizeof(double)*(N));
        post_pre = (double*) malloc(sizeof(double)*(N));
        
        //allocate space for final scores
        prior = (double*) malloc(sizeof(double)*(N));
        
        cout<<"read network done"<<endl;
        
    }

	void read_prior(){

		//initialize priors as theta_unl
		vertex node;
		for (node = 0; node < N; node++) {
			prior[node] = 0;
		}

		if (prior_file != "") {

			ifstream in(prior_file, ifstream::in);
			assert(in);

			string line;

			double score;
			while (getline(in, line) != NULL){
				node = (vertex)atol(strtok((char *)line.c_str(), " \n\t\r"));
				score = (double)atof(strtok(NULL, " \n\t\r"));
				prior[node] = score;
			}
			in.close();
		}

		if (train_set_file != ""){

			ifstream in(train_set_file, ifstream::in);
			assert(in);

			string line;

			//reading labeled sybil nodes.
			getline(in, line);
			getline(in, line);

			istringstream neg_train_str(line);
			vertex sub;
			while (neg_train_str){
				neg_train_str >> sub;
				neg_train_set.insert(sub);
			}

			set<vertex>::iterator iter;
			for (iter = neg_train_set.begin(); iter != neg_train_set.end(); iter++) {
				prior[*iter] = 1.0;
			}

			in.close();

		}

	}

    void write_posterior(){

        ofstream out(post_file, ofstream::out);

        for (vertex i = 0; i < N; i++) {
            out << i << " " << setprecision(10) << 1-post[i] << endl;
        }

        out.close();

    }

    static void * CIA_thread(void *arg_pointer){

        Data * pointer = ((CIA_arg *)arg_pointer)->data_pointer;
        int current_thread = ((CIA_arg *)arg_pointer)->current_thread;

        int num_nodes = ceil(((float)pointer->N) / pointer->num_threads);
        int start = current_thread * num_nodes;
        int end = current_thread * num_nodes + num_nodes;


        if (end > pointer->N) {
            end = pointer->N;
        }

        //cout << "begins " << current_thread << " " << start << " " << end << endl;

        vertex node;
        double message;

        //list<pair<vertex, double> >::iterator nei_iter;
        
        //pair<long, long> s;
        //stringstream convert1, convert2;
        
        set<pair<vertex, double> >::iterator iter;
        double nei_weight;
        double p;
        
        unordered_map<vertex,set<pair<vertex, double> > >::iterator iter_node, find_nei_wei;

        for (vertex index = start; index < end; index++) {
            node = pointer->ordering_array[index];
          
            message = 0;
            
			if (pointer->network_map.find(node) != pointer->network_map.end()){
				for (iter = pointer->network_map[node].begin(); iter != pointer->network_map[node].end(); iter++){
                    
                 
                    double sum_wei = 0;
					find_nei_wei = pointer->network_map.find((*iter).first);
					if (find_nei_wei != pointer->network_map.end()){
                        sum_wei = find_nei_wei->second.size();
                    }else{
                        cout<<"network out error at node "<<(*iter).first<<endl;
                    }
                
                    if(sum_wei == 0){
                        message += 0;
                    }else{
                        message += pointer -> post_pre[(*iter).first] * (*iter).second / sum_wei;
                    }
                }
            }
            
            pointer->post[node] = (1 - pointer->alpha) * message + pointer->alpha * pointer->prior[node];

        }

    }

 
    void CIA(){
        //use iterative method to propagate messages
        //flooding method's performance is not as good as iterative method

        ordering_array = (int*) malloc(sizeof(int) * (N) );

        //initialize posts
        memcpy(post, prior, sizeof(double) * (N));

        //random ordering
        vector<vertex> ordering;
        for (vertex i = 0; i < N; i++) {
            ordering.push_back(i);
        }

        int iter = 1;
        vector<vertex>::iterator iter_order;

        CIA_arg * arg_pointer;
        int current_thread;
        pthread_t thread;
        vector<pthread_t> threads;

        vertex i = 0;

        max_iter = (int)log(N);
        
        while (iter <= max_iter) {

            threads.clear();

            memcpy(post_pre, post, sizeof(double) * (N));

            random_shuffle(ordering.begin(), ordering.end(), p_myrandom);
            for (iter_order = ordering.begin(), i = 0; iter_order != ordering.end(); iter_order++, i++){
                ordering_array[i] = *iter_order;
            }

            for (current_thread = 0; current_thread < num_threads; current_thread++) {
                arg_pointer = new CIA_arg();
                arg_pointer->data_pointer = this;
                arg_pointer->current_thread = current_thread;

                pthread_create(&thread, NULL, CIA_thread, (void*)arg_pointer);
                threads.push_back(thread);
            }

            for (i = 0; i < threads.size(); i++) {
                pthread_join(threads[i], NULL);
                //cout << i << endl;
            }

          /*  //output result of this iteration
            stringstream post_file_iter;
            post_file_iter << post_file << iter << ".txt";
            const char *addr = post_file_iter.str().c_str();
            ofstream out(addr, ofstream::out);
            
            for (vertex i = 0; i < N; i++) {
                out << i << " " << setprecision(10) << post[i] << endl;
            }
            
            out.close(); */
            
            iter += 1;

        }
    }

    void parse_par(int argc, char **argv){

        //default setting
        network_file = "";
        alpha = 0.15;

        max_iter = 10;

        theta_pos = 0.9;
        theta_neg = 0.1;
        theta_unl = 0;

        train_set_file = "";
        post_file = "";
        prior_file = "";

        num_threads = 1;

		//by default, weighted graph
		weighted_graph = 1;
		weight = 0.9;


        int i = 1;

        while (i < argc) {

            if (strcmp(argv[i],"-graphfile") == 0){
                network_file = argv[i + 1];
            }
            else if (strcmp(argv[i],"-trainfile") == 0){
                train_set_file = argv[i + 1];
            }
            else if (strcmp(argv[i],"-priorfile") == 0){
                prior_file = argv[i + 1];
            }
            else if (strcmp(argv[i], "-postfile") == 0){
                post_file = argv[i + 1];
            }
            else if (strcmp(argv[i],"-mIter") == 0){
                max_iter = atoi(argv[i + 1]);
            }
			else if (strcmp(argv[i],"-alpha") ==0){
          	alpha = atof(argv[i + 1]);
			}
            else if (strcmp(argv[i],"-tp") == 0){
                theta_pos = atof(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-tn") == 0){
                theta_neg = atof(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-tu") == 0){
                theta_unl = atof(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-nt") == 0){
                num_threads = atoi(argv[i + 1]);
            }
            else{
                cout << "undefined inputs: " << argv[i] <<endl;
                exit(0);
            }

            i += 2;
        }

    }


};


int main (int argc, char **argv)
{

    srand ( time(NULL) );

    Data data;
    cout<<1<<endl;
    data.parse_par(argc, argv);
    cout<<2<<endl;
    data.read_network();
    cout<<3<<endl;
    data.read_prior();
    cout<<4<<endl;
    data.CIA();
	cout << 5 << endl;
    data.write_posterior();
  
    return 0;
}
