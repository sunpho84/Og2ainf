#ifndef INPUT_HPP
#define INPUT_HPP

#include "global.hpp"
#include "aliases.hpp"

enum ERR_t{NO_FAIL,FAILED_READ,FAILED_CONVERSION,FAILED_OPEN,MISPLACED_TK,UNINITIALIZED_PAR};

enum TK_glb_t{FEOF_GLB_TK,VALUE_GLB_TK, L_TK,APBC_TK,ACTION_TK,BETA_TK,CSW_TK,ALPHA_TK,R_TK,AL_TK};

const char L_tag[]="L";
const char APBC_tag[]="APBC";
const char action_tag[]="Action";
const char beta_tag[]="Beta";
const char csw_tag[]="Csw";
const char alpha_tag[]="Alpha";
const char r_tag[]="R";
const char al_tag[]="al";


// parse the value string
template <class T>
void get_value(FILE *fin,T &ret,const char *t);

// check
void check_str_par(const string str,const char *name);
void check_int_par(const int val,const char *name);
void check_double_par(const double val,const char *name);

// read global input file
void read_input_glb(const char input_path[]);

// read input file relative to single ensembles
// void read_input(const string &input, const string &name);

#endif
