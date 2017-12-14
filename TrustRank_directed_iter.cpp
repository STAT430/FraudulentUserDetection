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
#include <time.h>

#define GCC_VERSION (__GNUC__ * 10000 \
+ __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

#if GCC_VERSION >= 40300
#include <tr1/unordered_map>
using namespace std::tr1;
#define hash_map unordered_map

#else
#include <unordered_map>
#endif
using namespace std;

/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */

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
    unordered_map<vertex,set<pair<vertex, double> > > network_in, network_out, network_map;

    //1 : weighted, 0: unweighted
    int weighted_graph;

    //priors
    //probability of being benign for each node
    double *prior;

    //network file
    char *network_file;

    char *post_file;
    char *prior_file;


    //during computation, post is treated as the normalized multiplication of incoming messages of each node.
    //
    double *post;
    double *post_pre;

    //to support different input format
    char *train_set_file;


    //parameters of PMR model when input is a set of labeled nodes, instead of priors
    double theta_pos;
    double theta_neg;
    double theta_unl;

    //this weight is used when the input is an unweighted graph
    double weight;

    //for sybilrank
    double alpha;

    double max_iter;

    //-1: no boosting
    //0 : only benign labels are given. Sample the same number of nodes as Sybils
    //1 : only Sybils labels are givem. Sample the same number of nodes as benign labels
    //2 : assume no labels are given. Sample num_benign as benign labels. Sample num_sybil nodes as Sybil labels.
    int boosting_type;
    int num_benign;
    int num_sybil;

    int boosting_iter;

    int* ordering_array;

    int sampled_sybil_nodes;


    int num_threads;

    double polar;
    double homophily;

    double estimated_polar;

    string aggregator;


    class TrustRank_arg{
    public:
        Data * data_pointer;
        int current_thread;
        TrustRank_arg(){}
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
        pair <vertex, double> nei (node1, w);
        network_in[node2].insert(nei);
        network_out[node1].insert(make_pair(node2, w));
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
        vertex node1,node2,max_id=0;
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
            add_edge(node1, node2, w);
            if(max_id < node1){
                max_id = node1;
            }
            if(max_id < node2){
                max_id = node2;
            }
        }

        //number of nodes in the graph
//        N = network_map.size();
        N = max_id + 1;
//        cout<<N<<" "<<network_map.size()<<endl;

        //allocate space for final scores
        post = (double*) malloc(sizeof(double)*(N));
        post_pre = (double*) malloc(sizeof(double)*(N));

        //allocate space for final scores
        prior = (double*) malloc(sizeof(double)*(N));

        in.close();
        
        cout<<"read network done"<<endl;

    }


    void read_prior(){


        //initialize priors as theta_unl
        vertex node;
        
        if (train_set_file == "" && boosting_type == -1){
            for (node = 0; node < N; node++) {
                prior[node] = 0.5;
            }
        }else{
            for (node = 0; node < N; node++) {
                prior[node] = theta_unl;
            }
            cout<<"initialize prior done"<<endl;
        }

        if (prior_file != "") {


            ifstream in(prior_file,ifstream::in);
            assert(in);

            string line;

            double score;
            while(getline(in,line)!=NULL){
                node=(vertex)atol(strtok((char *)line.c_str()," \n\t\r"));
                score=(double)atof(strtok(NULL," \n\t\r"));
                prior[node] = score;
            }
            in.close();
        }

        if (train_set_file != "" && boosting_type == -1){

            ifstream in(train_set_file,ifstream::in);
            assert(in);


            string line;

            //reading labeled benign nodes.
            getline(in,line);
            istringstream pos_train_str(line);
            vertex sub;
            while (pos_train_str){
                pos_train_str >> sub;
                prior[sub] = theta_pos;
                //cout << sub << " " << prior[sub] << endl;
            }

            getline(in,line);
            istringstream neg_train_str(line);
            // while (neg_train_str){
            //     neg_train_str >> sub;
            //     prior[sub] = theta_neg;
            //     //cout << sub << " " << prior[sub] << endl;
            // }

            in.close();

        }
        
        cout<<"read prior done"<<endl;

    }


    void write_posterior(){

        ofstream out(post_file, ofstream::out);

        for (vertex i = 0; i < N; i++) {
            out << i << " " << setprecision(10) << post[i] << endl;
        }

        out.close();

    }

    //output peak memory use and total runtime
    // void write_memtime(long mem, double time){
    //
    //     ofstream out(post_file, ofstream::out);
    //     out << "Memory\t" << mem << endl;
    //     out << "Time\t" << time << endl;
    //
    //     out.close();
    //     cout<< "Memory\t" << mem << endl;
    //     cout<< "Time\t" << time << endl;
    //
    // }

    static void * TrustRank_thread(void *arg_pointer){

        Data * pointer = ((TrustRank_arg *)arg_pointer)->data_pointer;
        int current_thread = ((TrustRank_arg *)arg_pointer)->current_thread;

        int num_nodes = ceil(((float)pointer->N) / pointer->num_threads);
        int start = current_thread * num_nodes;
        int end = current_thread * num_nodes + num_nodes;


        if (end > pointer->N) {
            end = pointer->N;
        }

        //cout << "begins " << current_thread << " " << start << " " << end << endl;


        vertex node;
        double message;

        set<pair<vertex, double> >::iterator nei_iter;

        double value[2], sum;

//## Code below edited by Le Zhang on Feb 16, 2016 ##
//## lezhang@iastate.edu

        for (vertex index = start; index < end; index++) {
            node = pointer->ordering_array[index];

            //update the the post for node
            message = 0;
            if(pointer->network_in.find(node) != pointer->network_in.end()){
                for (nei_iter = pointer->network_in[node].begin(); nei_iter != pointer->network_in[node].end(); nei_iter++) {
                    
                    double sum_wei = 0;
                    if(pointer->network_out.find((*nei_iter).first) != pointer->network_out.end()){
                        sum_wei = pointer->network_out[(*nei_iter).first].size();
                    }else{
                        cout<<"network out error at node "<<(*nei_iter).first<<endl;
                    }
                    //use edge weight but will be extremely slow
                    //                if(pointer->network_out[(*nei_iter).first].size()!=0){
                    //                    for(set<pair<vertex, double> >::iterator wei_iter = pointer->network_out[(*nei_iter).first].begin(); wei_iter!=pointer->network_out[(*nei_iter).first].end(); wei_iter++){
                    //                        sum_wei += (*wei_iter).second;
                    //                    }
                    //                }
                    if(sum_wei == 0){
                        message += 0;
                    }else{
                        message += pointer -> post_pre[(*nei_iter).first] * (*nei_iter).second / sum_wei;
                    }
                    
                    
                    // * (*nei_iter).second;
                    // cout<<message<<endl;
                    
                }
            }
            
            pointer->post[node] = (1 - pointer->alpha) * message + pointer->alpha * pointer->prior[node];

        }

        //cout << "end " << current_thread << " " << start << " " << end << endl;

    }


    void TrustRank(){
        //use iterative method to propagate messages
        //flooding method's performance is not as good as iterative method

        ordering_array = (int*) malloc(sizeof(int) * N);

        //initialize posts
        memcpy(post, prior, sizeof(double) * N);

        //random ordering
        vector<vertex> ordering;
        for (vertex i = 0; i < N; i++) {
            ordering.push_back(i);
        }

        int iter = 1;
        vector<vertex>::iterator iter_order;

        TrustRank_arg * arg_pointer;
        int current_thread;
        pthread_t thread;
        vector<pthread_t> threads;

        vertex i = 0;


        while (iter <= max_iter) {

            threads.clear();

            memcpy(post_pre, post, sizeof(double) * N);

            random_shuffle(ordering.begin(), ordering.end(), p_myrandom);
            for (iter_order = ordering.begin(), i = 0; iter_order != ordering.end(); iter_order++, i++){
                ordering_array[i] = *iter_order;
            }

            for (current_thread = 0; current_thread < num_threads; current_thread++) {
                arg_pointer = new TrustRank_arg();
                arg_pointer->data_pointer = this;
                arg_pointer->current_thread = current_thread;

                pthread_create(&thread, NULL, TrustRank_thread, (void*)arg_pointer);
                threads.push_back(thread);
            }

            for (i = 0; i < threads.size(); i++) {
                pthread_join(threads[i], NULL);
                //cout << i << endl;
            }

			//output result of this iteration

			//post_file_iter = post_file + to_string(iter) + ".txt";
			stringstream post_file_iter;
			post_file_iter << post_file << iter << ".txt";
			const char *addr = post_file_iter.str().c_str();
			ofstream out(addr, ofstream::out);

			for (vertex i = 0; i < N; i++) {
				out << i << " " << setprecision(10) << post[i] << endl;
			}

			out.close();

			iter += 1;

			cout << "write graph at iter " << iter << "done." << endl;

        }

        // Normalization
//        for(vertex v=0; v<N;v++){
//            if(network_in[v].size() != 0){
//                post[v] /= network_in[v].size();
//            }
//        }

    }

    void parse_par(int argc, char **argv){

        //default setting
        network_file = "";
        boosting_type = -1;
        boosting_iter = 5;

        num_benign = 10;
        num_sybil = 10;


        max_iter = 10;


        theta_pos = 1.0;
        theta_neg = 0.1;
        theta_unl = 0;

        weight = 1.0;

        alpha = 0.15;//for sybilrank

        train_set_file = "";
        post_file = "";
        prior_file = "";

        //by default, weighted graph
        weighted_graph = 1;

        sampled_sybil_nodes = 1;

        num_threads = 1;

        aggregator = "hea";

        estimated_polar = -1;

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
            else if (strcmp(argv[i],"-bt") == 0){
                boosting_type = atoi(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-bi") == 0){
                boosting_iter = atoi(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-mIter") == 0){
                max_iter = atoi(argv[i + 1]);
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
            else if (strcmp(argv[i],"-wei") == 0){
                weight = atof(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-wg") == 0){
                weighted_graph = atoi(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-ss") == 0){
                sampled_sybil_nodes = atoi(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-nt") == 0){
                num_threads = atoi(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-nb") == 0){
                num_benign = atoi(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-ns") == 0){
                num_sybil = atoi(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-agg") == 0){
                aggregator = argv[i + 1];
            }
            else if (strcmp(argv[i],"-ep") == 0){
                estimated_polar = atof(argv[i + 1]);
            }
            else if (strcmp(argv[i],"-alph") == 0){
                alpha = atof(argv[i + 1]);
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
    clock_t start, end;
    start = clock() ;

    srand ( time(NULL) );

    Data data;
	cout << 1 << endl;
    data.parse_par(argc, argv);
	cout << 2 << endl;
    data.read_network();
	cout << 3 << endl;
    data.read_prior();
	cout << 4 << endl;
	data.TrustRank();
    // data.write_posterior();
   
    end = clock() ;
//    cout<<endl<<"Total time taken: "<<(double)(end-start)/CLOCKS_PER_SEC*1000<<" ms"<<endl;
//    size_t peakSize = getPeakRSS( );
//    cout<<peakSize<<endl;
    // data.write_memtime((long)peakSize, (double)(end-start)/CLOCKS_PER_SEC*1000);

    return 0;
}
