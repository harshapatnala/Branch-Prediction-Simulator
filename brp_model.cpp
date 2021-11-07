#include <iostream>
#include <cmath>
#include <cstdio>
#include "brp_model.h"
using namespace std;

bimodal_bp::bimodal_bp(unsigned int num_pc_index_bits) { //Constructor
    index_bits = num_pc_index_bits;
    mispredictions = predictions = index_mask = 0;
    misprediction_rate = 0;
    pc_count = 0;
    num_sets = (unsigned int) pow(2, index_bits);
    counter_array = new int[num_sets];
    for(unsigned int i=0; i < num_sets; i++) {
        counter_array[i] = 2;
    }
    for(unsigned int i=0; i < index_bits; i++) {
        index_mask <<= 1;
        index_mask +=1;
    }
}

bimodal_bp::~bimodal_bp() { //Destructor
    delete counter_array;
}
char bimodal_bp::predicted_outcome(unsigned long int program_counter, unsigned int & set_number) {
    program_counter >>= 2;
    set_number = program_counter & index_mask;
    int counter = counter_array[set_number];
    if(counter >=2) return 't';
    else return 'n';
}

void bimodal_bp::update_counter(unsigned int set_number, char predict_outcome, char actual_outcome) {
    int counter = counter_array[set_number];
    if(predict_outcome=='t' && actual_outcome == 't') {
        predictions++;
        if(counter <3) counter_array[set_number]++;
    }
    else if(predict_outcome=='n' && actual_outcome =='n') {
        predictions++;
        if(counter >0) counter_array[set_number]--;
    }

    else if(predict_outcome =='t' && actual_outcome == 'n') {
        mispredictions++;
        counter_array[set_number]--;
    }
    else if(predict_outcome =='n' && actual_outcome == 't') {
        mispredictions++;
        counter_array[set_number]++;
    }
}
void bimodal_bp::predict_outcome(unsigned long int program_counter, char actual_outcome) {
    pc_count++;
    unsigned int set_number;
    char p_outcome = predicted_outcome(program_counter, set_number);
    update_counter(set_number, p_outcome, actual_outcome);
}

void bimodal_bp::print_contents() {
    printf("FINAL BIMODAL CONTENTS\n");
    for(unsigned int i=0; i < num_sets; i++) {
        printf("%d\t%d\n", i, counter_array[i]);
    }
}
void bimodal_bp::print_stats() {
    misprediction_rate = (double) mispredictions/ pc_count *100;
    printf("OUTPUT\n");
    printf("number of predictions:     %d\n", pc_count);
    printf("number of mispredictions:  %d\n", mispredictions);
    printf("misprediction rate:        %.2f%%\n", misprediction_rate);
}


gshare_bp::gshare_bp(unsigned int m1, unsigned int m2) { // Constructor
    pc_index_bits = m1;
    gbhr_bits = m2;
    pc_index_mask = gbhr_mask = mispredictions = predictions = pc_count = gbhr= pc_lower_offset = pc_upper_mask =pc_lower_mask=0;
    misprediction_rate = 0;
    gbhr_taken_mask = 1;
    gbhr_taken_mask <<= (int)(gbhr_bits-1);

    num_sets = (unsigned int) pow(2, pc_index_bits);
    counter_array = new int[num_sets];
    for(unsigned int i=0; i < num_sets; i++) {
        counter_array[i] = 2;
    }
    for(unsigned int i=0; i < gbhr_bits; i++) {
        gbhr_mask <<=1;
        gbhr_mask += 1;
    }
    for(unsigned int i=0; i < pc_index_bits; i++) {
        pc_index_mask <<= 1;
        pc_index_mask += 1;
    }

    for(unsigned int i=0; i < (pc_index_bits - gbhr_bits); i++) {
        pc_lower_mask <<= 1;
        pc_lower_mask +=1;
    }
}

gshare_bp::~gshare_bp() { //Destructor
    delete counter_array;
}
char gshare_bp::predicted_outcome(unsigned long int program_counter, unsigned int & set_number) {
    program_counter >>= 2;
    program_counter = program_counter & pc_index_mask;
    pc_lower_offset = program_counter & pc_lower_mask;
    pc_upper_mask = program_counter >> (pc_index_bits - gbhr_bits);
    gbhr &= gbhr_mask;
    pc_upper_mask = pc_upper_mask ^ gbhr;
    pc_upper_mask <<= (pc_index_bits - gbhr_bits);
    set_number = pc_upper_mask | pc_lower_offset;
    unsigned int count = counter_array[set_number];
    if(count >=2) return't';
    else return'n';
}

void gshare_bp::update_counter(unsigned int set_number, char p_outcome, char actual_outcome) {
    unsigned int count = counter_array[set_number];
    if(p_outcome == 't' && actual_outcome == 't') {
        predictions++;
        if(count <3) counter_array[set_number]++;
    }
    else if(p_outcome == 'n' && actual_outcome == 'n') {
        predictions++;
        if(count >0) counter_array[set_number]--;
    }
    else if(p_outcome == 't' && actual_outcome == 'n') {
        mispredictions++;
        counter_array[set_number]--;
    }
    else if(p_outcome =='n' && actual_outcome == 't') {
        mispredictions++;
        counter_array[set_number]++;
    }
}

void gshare_bp::update_gbhr(char actual_outcome) {
    gbhr >>=1;
    if(actual_outcome=='t') gbhr |= gbhr_taken_mask;
}
void gshare_bp::predict_outcome(unsigned long int program_counter, char actual_outcome) {

    pc_count++;
    unsigned int set_number;
    char p_outcome = predicted_outcome(program_counter, set_number);
    update_counter(set_number, p_outcome, actual_outcome);
    update_gbhr(actual_outcome);
}

void gshare_bp::print_contents() {
    printf("FINAL GSHARE CONTENTS\n");
    for(unsigned int i=0; i < num_sets; i++) {
        printf("%d\t%d\n", i, counter_array[i]);
    }
}
void gshare_bp::print_stats() {
    misprediction_rate = (double)mispredictions/ pc_count * 100;
    printf("OUTPUT\n");
    printf("number of predictions:     %d\n", pc_count);
    printf("number of mispredictions:  %d\n", mispredictions);
    printf("misprediction rate:        %.2f%%\n", misprediction_rate);
}


hybrid_bp::hybrid_bp(unsigned int c_bits, unsigned int g_pc_bits, unsigned int g_gb_bits, unsigned int b_bits) {
    chooser_bits = c_bits;
    chooser_num_sets = (unsigned int) pow(2, c_bits);
    pc_count = pc_index_mask = mispredictions = 0;
    misprediction_rate = 0;
    chooser_array = new int[chooser_num_sets];
    gshare_predictor = new gshare_bp(g_pc_bits, g_gb_bits);
    bimodal_predictor = new bimodal_bp(b_bits);
    for(unsigned int i=0; i < chooser_num_sets; i++){
        chooser_array[i] = 1;
    }
    for(unsigned int i=0; i < chooser_bits; i++) {
        pc_index_mask <<= 1;
        pc_index_mask += 1;
    }
}

void hybrid_bp::predict_outcome(unsigned long int program_counter, char actual_outcome) {
    pc_count++;
    unsigned int bimodal_set_number, gshare_set_number;
    char bimodal_predicted_outcome = bimodal_predictor->predicted_outcome(program_counter, bimodal_set_number);
    char gshare_predicted_outcome = gshare_predictor->predicted_outcome(program_counter, gshare_set_number);
    program_counter >>= 2;
    unsigned int chooser_set_number = program_counter & pc_index_mask;
    unsigned int chooser_count = chooser_array[chooser_set_number];
    char chooser_predicted_outcome;

    if(chooser_count >=2) {
        chooser_predicted_outcome = gshare_predicted_outcome;
        gshare_predictor->update_counter(gshare_set_number, chooser_predicted_outcome, actual_outcome);
    }
    else {
        chooser_predicted_outcome = bimodal_predicted_outcome;
        bimodal_predictor->update_counter(bimodal_set_number, chooser_predicted_outcome, actual_outcome);
    }
    if(chooser_predicted_outcome != actual_outcome) mispredictions++;

    if((bimodal_predicted_outcome==actual_outcome) && (gshare_predicted_outcome !=actual_outcome)) {
        if(chooser_count >0) chooser_array[chooser_set_number]--;
    }
    else if((gshare_predicted_outcome == actual_outcome) && (bimodal_predicted_outcome != actual_outcome)) {
        if(chooser_count < 3) chooser_array[chooser_set_number]++;
    }
    gshare_predictor->update_gbhr(actual_outcome);
}
void hybrid_bp::print_contents() {
    printf("FINAL CHOOSER CONTENTS\n");
    for(unsigned int i=0; i < chooser_num_sets; i++) {
        printf("%d\t%d\n", i, chooser_array[i]);
    }
    gshare_predictor->print_contents();
    bimodal_predictor->print_contents();
}

void hybrid_bp::print_stats() {
    misprediction_rate = (double) mispredictions / pc_count *100;
    printf("OUTPUT\n");
    printf("number of predictions:     %d\n", pc_count);
    printf("number of mispredictions:  %d\n", mispredictions);
    printf("misprediction rate:        %.2f%%\n", misprediction_rate);
}

hybrid_bp::~hybrid_bp() {
    delete chooser_array;
    delete bimodal_predictor;
    delete gshare_predictor;
}