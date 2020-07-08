#include <iostream>
#include "read_input.hpp"
#include "aliases.hpp"
#include "global.hpp"
// #include "operations.hpp"
// #include "read.hpp"

#define DEFAULT_STR_VAL "null"
#define DEFAULT_INT_VAL -1
#define DEFAULT_DOUBLE_VAL 1.2345

using namespace std;


char tok[128];

TK_glb_t get_TK_glb(FILE *fin)
{
    //read a token
    int rc=fscanf(fin,"%s",tok);
    if(rc!=1)
    {
        if(feof(fin)) return FEOF_GLB_TK;
        else
        {
            fprintf(stderr,"Getting %d while reading token\n",rc);
            exit(FAILED_READ);
        }
    }

    //parse the token
    if(strcasecmp(tok,L_tag)==0) return L_TK;
    if(strcasecmp(tok,APBC_tag)==0) return APBC_TK;
    if(strcasecmp(tok,action_tag)==0) return ACTION_TK;
    if(strcasecmp(tok,beta_tag)==0) return BETA_TK;
    if(strcasecmp(tok,csw_tag)==0) return CSW_TK;
    if(strcasecmp(tok,alpha_tag)==0) return ALPHA_TK;
    if(strcasecmp(tok,r_tag)==0) return R_TK;
    if(strcasecmp(tok,al_tag)==0) return AL_TK;

    return VALUE_GLB_TK;
}


//parse the value string (glb input)
template <class T>
void _get_value_glb(FILE *fin,T &ret,const char *t)
{
    TK_glb_t tk=get_TK_glb(fin);
    if(tk!=VALUE_GLB_TK)
    {
        fprintf(stderr,"Getting token %s in the wrong place\n",tok);
        exit(MISPLACED_TK);
    }

    int rc=sscanf(tok,t,&ret);

    if(rc!=1)
    {
        fprintf(stderr,"Converting %s to %s failed\n",tok,t);
        exit(FAILED_CONVERSION);
    }
}

void get_value_glb(FILE *fin,double &out)
{
    return _get_value_glb(fin,out,"%lg");
}

void get_value_glb(FILE *fin,int &out)
{
    return _get_value_glb(fin,out,"%d");
}

void get_value_glb(FILE *fin,string &out)
{
    char temp[1024];
    _get_value_glb(fin,temp,"%s");
    out=string(temp);

}


//check
void check_str_par(const string str,const char *name)
{
    if(str.compare(DEFAULT_STR_VAL)==0)
    {
        fprintf(stderr,"%s not initialized\n",name);
        exit(UNINITIALIZED_PAR);
    }
}

void check_int_par(const int val,const char *name)
{
    if(val==DEFAULT_INT_VAL)
    {
        fprintf(stderr,"%s not initialized\n",name);
        exit(UNINITIALIZED_PAR);
    }
}

void check_double_par(const double val,const char *name)
{
    if(val==DEFAULT_DOUBLE_VAL)
    {
        fprintf(stderr,"%s not initialized\n",name);
        exit(UNINITIALIZED_PAR);
    }
}

// reads the input file
void read_input_glb(const char path[])
{
    FILE *fin=fopen(path,"r");
    if(not fin)
    {
        fprintf(stderr,"Failed to open \"%s\"\n",path);
        exit(FAILED_OPEN);
    }

    L      = DEFAULT_INT_VAL;
    APBC   = DEFAULT_INT_VAL;
    action = DEFAULT_INT_VAL;
    beta   = DEFAULT_DOUBLE_VAL;
    csw    = DEFAULT_DOUBLE_VAL;
    alpha  = DEFAULT_DOUBLE_VAL;
    r      = DEFAULT_DOUBLE_VAL;
    al     = DEFAULT_DOUBLE_VAL;

    while(not feof(fin))
    {
        TK_glb_t tk=get_TK_glb(fin);
        switch(tk)
        {
            case VALUE_GLB_TK:
                fprintf(stderr,"Invalid token %s found\n",tok);
                exit(1);
                break;
            case L_TK:
                get_value_glb(fin,L);
                break;
            case APBC_TK:
                get_value_glb(fin,APBC);
                break;
            case ACTION_TK:
                get_value_glb(fin,action);
                break;
            case BETA_TK:
                get_value_glb(fin,beta);
                break;
            case CSW_TK:
                get_value_glb(fin,csw);
                break;
            case ALPHA_TK:
                get_value_glb(fin,alpha);
                break;
            case R_TK:
                get_value_glb(fin,r);
                break;
            case AL_TK:
                get_value_glb(fin,al);
                break;

            case FEOF_GLB_TK:
                break;
        }
    }

    //check initialization
    check_int_par(L, L_tag);
    check_int_par(APBC, APBC_tag);
    check_int_par(action,action_tag);
    check_double_par(beta,beta_tag);
    check_double_par(csw,csw_tag);
    check_double_par(alpha,alpha_tag);
    check_double_par(r,r_tag);
    check_double_par(al,al_tag);

    fclose(fin);

    //print input parameters
    printf("*------------------------------------------------------*\n");
    printf("|                Global configuration                  |\n");
    printf("*------------------------------------------------------*\n\n");

    T = 2*L;
    dim << L,T;
    V << L,L,L,T;
    printf(" Volume = %d x %d\n",L,T);

    printf("    with BC: ");
    if(APBC) printf("AntiPeriodic\n\n");
    else     printf("Periodic\n\n");

    string action_name;
    if(action==1)     action_name =  "Wilson";
    else if(action==2)action_name =  "Tree Level Symanzik";
    else if(action==9)action_name =  "Iwasaki";
    else {cout<<"Action not available"<<endl<<endl; exit(0);}

    c1 = (action==1 ? 0.0 : (action==2 ? -1.0/12.0 : (action==9 ? -0.331 : (1.0/0.0))));

    printf("Action: %s  (c1 = %.3lf)",action_name.c_str(),c1);


    printf("\n\n");
}
