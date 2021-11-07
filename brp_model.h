
#ifndef PROJ2_BRANCH_SIM_BRP_MODEL_H
#define PROJ2_BRANCH_SIM_BRP_MODEL_H


class bimodal_bp {
protected:
    int* counter_array;
    unsigned int index_bits, num_sets, index_mask, mispredictions, predictions, pc_count;
    double misprediction_rate;

public:
    bimodal_bp(unsigned int);
    ~bimodal_bp();

    char predicted_outcome(unsigned long int, unsigned int&);
    void update_counter(unsigned int, char, char);
    void predict_outcome(unsigned long int, char);
    void print_stats();
    void print_contents();
};

class gshare_bp {
protected:
    int* counter_array;
    unsigned int num_sets, pc_index_bits, gbhr_bits, pc_index_mask, pc_upper_mask, gbhr, gbhr_mask,
    pc_lower_offset, pc_lower_mask, mispredictions, predictions, pc_count;
    int gbhr_taken_mask;
    double misprediction_rate;

public:
    gshare_bp(unsigned int, unsigned int);
    ~gshare_bp();

    char predicted_outcome(unsigned long int, unsigned int&);
    void update_counter(unsigned int, char, char);
    void update_gbhr(char);
    void predict_outcome(unsigned long int, char);
    void print_contents();
    void print_stats();
};


class hybrid_bp {
protected:
    int* chooser_array;
    unsigned int chooser_bits, chooser_num_sets, pc_count, pc_index_mask, mispredictions;
    double misprediction_rate;
    bimodal_bp* bimodal_predictor;
    gshare_bp* gshare_predictor;
public:
    hybrid_bp(unsigned int, unsigned int, unsigned int, unsigned int);
    ~hybrid_bp();
    void predict_outcome(unsigned long int, char);
    void print_stats();
    void print_contents();
};


#endif //PROJ2_BRANCH_SIM_BRP_MODEL_H
