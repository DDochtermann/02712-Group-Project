//
//  Memory and changing population sizes -- here with waiting time till new mutation introduces the new allele
//  Basically now have a small ate of mutation possible from A to a (but not reverse mutation)

//
//  Copyright © 2016 Oana. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//#include <malloc.h>

/*
 *  Special global variable, which stores the kind of uniform random number generator we use from gsl
 */
const gsl_rng *gBaseRand;

/*
 *  Simulation parameters.
 */
#define NRRUNS 10000    // Total number of runs in this batch.
#define G 100000000    // Maximum number of generations.

#define MUTANT     1
#define WILDTYPE     0

#define DEATH_RATE_W     1 //intrinsic to the genotype, or phenotype
#define DEATH_RATE_M     1



/************/
/*  Main    */
/************/



int main(int argc, const char * argv[]) {
    clock_t begin = clock(), end; //begin clock
    double time_spent;
    
    int K =5000; //Assign maximum carrying capacity //atoi(argv[1]);
    int N_0 =1000; //Assign initial population size to wildtype population //atoi(argv[1]);
    
    
    // double BIRTH_RATE_W_E1=atof(argv[4]);
    double BIRTH_RATE_W_E1=0.95;
    double variance_e = atof(argv[1]);
    
    double BIRTH_RATE_M1_E1=BIRTH_RATE_W_E1 + variance_e;
    // double BIRTH_RATE_M1_E1=atof(argv[1]);
    
    double BIRTH_RATE_M2_E1=BIRTH_RATE_W_E1 - variance_e;
    // double BIRTH_RATE_M2_E1=atof(argv[2]);
    
    
    double p=atof(argv[2]); //memory, eventually as argument
    double mutationAa=0.00001;
    
    int i=0, j=0, indiv=0;  // General loop counter.
    int run=0, counter=0;    // Loop counter for number of runs of the simulation in this batch.
    
    /*  record */
    
    int threshhold=0,
    extinctions=0;
    
    /*  record */
    double timetoextinctall[NRRUNS];
    double timetoextinctgenerations[NRRUNS];
    double timetoextinctgenerationsdeaths[NRRUNS];
    double timeofintroduction[NRRUNS];
    int popsizeatintroduction[NRRUNS];
    double sumtimetoextinct=0.0f,
    meantimetoextinct=0.0f,
    vartimetoextinct=0.0f;
    
    FILE  *outfile, *outfile2;
    
    struct timeval time; //using time to seed every run; see below for details (within the run)
    /* specifying to use Mersenne twister MT-19937 as the uniform PRNG */
    gBaseRand = gsl_rng_alloc(gsl_rng_mt19937);
    
    
    setbuf (stdout, NULL); /* Turn off buffering on stdout so         */
    /*     progress reports appear properly.   */
    
    for (run=0; run<NRRUNS; run++){
        //over RUNS HERE
        // printf("run %d\n", run);
        
        /*
         *  Set seed for this run.
         */
        gettimeofday(&time, NULL);
        unsigned   int randSeed = time.tv_sec*1000000+(time.tv_usec);  /* returns a non-negative integer */
        //printf ("  seed %lu \n", randSeed);
        gsl_rng_set (gBaseRand, randSeed);    /* seed the PRNG */
        
        
        /*
         *  Dynamic attributes.
         */
        int N = N_0; //Initialize the changing population size.
        //printf("initial N is %d\n", N);
        int W_Size = N_0; //Initialize variable for wildtype population size tracking
        int M1_Size = 0; //Initialize variable for mutant population size tracking
        int M2_Size = 0; //Initialize variable for mutant population size tracking
        int M_Size=M1_Size+M2_Size;
        int NR_Births=0, NR_Deaths=0, TOTAL_NR_Events=0;
        double firstIntroduction=0.0f;
        //int nrMutantsA=0; //nr of mutants at second locus for A
        //int nrMutantsM1=0;//nr of mutants at second locus for M1
        //int nrMutantsM2=0;//nr of mutants at second locus for M2
        
        
        
        // these new mutants basically accumulate with every birth, every single time A or M1 or M2 is born, there is a probability of mutation at this second locus
        
        //no introduction of a new mutant here -- allow mutaion between A and a and that will introduce it
        //int rnd;    //add one mutant: first select which a.
        //rnd = gsl_rng_uniform_int(gBaseRand,2); //put zero here if you want to test without a mutant introduction
        //rnd=1;
        //if(rnd==0){
        //    M1_Size++;
        //} else if(rnd==1){
        //    M2_Size++;
        //}
        //M_Size=M1_Size+M2_Size;
        //N=W_Size+M1_Size+M2_Size;
        
        
        //have a vector here the size of the carrying cappacity K
        // for the alleleic profiles
        
        
        int theEnv=1; //always start with environment 1 for now (bad for A -- drug)
        double B_rate_A=BIRTH_RATE_W_E1*W_Size; //*(1-N/K);
        double D_rate_A=DEATH_RATE_W*W_Size;
        double B_rate_M1=BIRTH_RATE_M1_E1*M1_Size; //*(1-N/K);
        double B_rate_M2=BIRTH_RATE_M2_E1*M2_Size; //*(1-N/K);
        double D_rate_M1=DEATH_RATE_M*M1_Size;
        double D_rate_M2=DEATH_RATE_M*M2_Size;
        double theRates[6]={B_rate_A, D_rate_A, B_rate_M1, B_rate_M2, D_rate_M1, D_rate_M2};
        double SUM_Rates=B_rate_A+D_rate_A+B_rate_M1+B_rate_M2+D_rate_M1+D_rate_M2;
        
        
        double timetoextinct=0.0f,
        event_time_total=0.0f,
        current_env_sum_time=0.0f,
        next_time=0.0f;
        
        counter=0;
        int havebreaked=0;
        
        while (N>0 ){
            counter++;
            if (N>=(K-1)){
                havebreaked=1;
                printf("breaking after 10000000 generations \n");
                printf("breaked W %d\n", W_Size);
                printf("breaked m1 %d\n", M1_Size);
                printf("breaked m2 %d\n", M2_Size);
                break;
            }
            
            SUM_Rates=B_rate_A+D_rate_A+B_rate_M1+B_rate_M2+D_rate_M1+D_rate_M2;
            next_time=gsl_ran_exponential (gBaseRand, 1.0/SUM_Rates); //Draw the waiting time until an event occurs from an exponential distribution.
            event_time_total=event_time_total+next_time; //Add the waiting time to the current (previous event) time
            current_env_sum_time=current_env_sum_time+next_time;
            
            
            //Sample event.
            double partialSum=0.0f;
            int theEventNr=0;
            
            double randomNumber2 = (gsl_rng_uniform_pos(gBaseRand) * SUM_Rates);
            
            theRates[0]=B_rate_A;
            theRates[1]=D_rate_A;
            theRates[2]=B_rate_M1;
            theRates[3]=B_rate_M2;
            theRates[4]=D_rate_M1;
            theRates[5]=D_rate_M2;
            
            int rate=0;
            while (randomNumber2 > partialSum)
            {
                partialSum += theRates[rate];
                rate++;
            }
            theEventNr=rate-1;
            TOTAL_NR_Events++;
            
            if (theEventNr==0 & W_Size>0){ //B_rate_A
                NR_Births++;
    
                    
                // here need to add a (constant) probability of mutation to a
                unsigned  int mutate=gsl_ran_binomial (gBaseRand, mutationAa, 1); //probability to genetically mutate
                if (mutate==1 ){ //1 = mutates;
                    if (firstIntroduction==0){
                        firstIntroduction=1;
                        timeofintroduction[run]=event_time_total;
                        printf("time %f\n", event_time_total);
                        popsizeatintroduction[run]=W_Size;
                        printf("W_Size %d\n", W_Size);
                    }
                    // M1_Size++;
                    float p_mutate = (float)(M2_Size+1)/((float)(M1_Size+1)+(float)(M2_Size+1));
                    unsigned  int which=gsl_ran_binomial (gBaseRand, p_mutate, 1);
                    if (which==0){
                        M1_Size++;
                    
                    } else {
                        M2_Size++;                        
                    }
                } else {
                    W_Size++; //else keep the A allele
                }
                
                N=W_Size+M1_Size+M2_Size;
                if (theEnv==1){
                    B_rate_A=BIRTH_RATE_W_E1*W_Size; //*(1-N/K);
                    D_rate_A=DEATH_RATE_W*W_Size;
                    B_rate_M1=BIRTH_RATE_M1_E1*M1_Size; //*(1-N/K);
                    B_rate_M2=BIRTH_RATE_M2_E1*M2_Size; //*(1-N/K);
                    D_rate_M1=DEATH_RATE_M*M1_Size;
                    D_rate_M2=DEATH_RATE_M*M2_Size;
                    SUM_Rates=B_rate_A+D_rate_A+B_rate_M1+B_rate_M2+D_rate_M1+D_rate_M2;
                }
    
                
            } else if (theEventNr==1 & W_Size>0){ //death of an A
                NR_Deaths++;
                W_Size--;
                N=W_Size+M1_Size+M2_Size;
                if (theEnv==1){
                    B_rate_A=BIRTH_RATE_W_E1*W_Size; //*(1-N/K);
                    D_rate_A=DEATH_RATE_W*W_Size;
                    B_rate_M1=BIRTH_RATE_M1_E1*M1_Size; //*(1-N/K);
                    B_rate_M2=BIRTH_RATE_M2_E1*M2_Size; //*(1-N/K);
                    D_rate_M1=DEATH_RATE_M*M1_Size;
                    D_rate_M2=DEATH_RATE_M*M2_Size;
                    SUM_Rates=B_rate_A+D_rate_A+B_rate_M1+B_rate_M2+D_rate_M1+D_rate_M2;
                }
                
                
                
            } else if (theEventNr==2 & M1_Size>0){ //Birth M1
                //printf("B M1:\n ");
                NR_Births++;

                unsigned  int tokeep=gsl_ran_binomial (gBaseRand, p, 1);
                if (tokeep==0){ //0 = redraw phenotype.
                    unsigned  int which=gsl_rng_uniform_int(gBaseRand,2);
                    if (which==0){
                        M1_Size++;
                    } else {
                        M2_Size++;
                    }
                } else {
                    M1_Size++;
                }
                
                M_Size=M1_Size+M2_Size;
                N=W_Size+M1_Size+M2_Size;
                
                if (theEnv==1){
                    B_rate_A=BIRTH_RATE_W_E1*W_Size; //*(1-N/K);
                    D_rate_A=DEATH_RATE_W*W_Size;
                    B_rate_M1=BIRTH_RATE_M1_E1*M1_Size; //*(1-N/K);
                    B_rate_M2=BIRTH_RATE_M2_E1*M2_Size; //*(1-N/K);
                    D_rate_M1=DEATH_RATE_M*M1_Size;
                    D_rate_M2=DEATH_RATE_M*M2_Size;
                    SUM_Rates=B_rate_A+D_rate_A+B_rate_M1+B_rate_M2+D_rate_M1+D_rate_M2;
                }
                
                
                
            } else if (theEventNr==3 & M2_Size>0){ //B M2
                NR_Births++;
                    
                unsigned  int tokeep=gsl_ran_binomial (gBaseRand, p, 1);
                if (tokeep==0){ //0 = redraw phenotype.
                    unsigned  int which=gsl_rng_uniform_int(gBaseRand,2);
                    if (which==0){
                        M1_Size++;
                    } else {
                        M2_Size++;
                    }
                } else {
                    M2_Size++;
                }
                
                M_Size=M1_Size+M2_Size;
                N=W_Size+M1_Size+M2_Size;
                
                if (theEnv==1){
                    B_rate_A=BIRTH_RATE_W_E1*W_Size; //*(1-N/K);
                    D_rate_A=DEATH_RATE_W*W_Size;
                    B_rate_M1=BIRTH_RATE_M1_E1*M1_Size; //*(1-N/K);
                    B_rate_M2=BIRTH_RATE_M2_E1*M2_Size; //*(1-N/K);
                    D_rate_M1=DEATH_RATE_M*M1_Size;
                    D_rate_M2=DEATH_RATE_M*M2_Size;
                    SUM_Rates=B_rate_A+D_rate_A+B_rate_M1+B_rate_M2+D_rate_M1+D_rate_M2;
                }
                
            } else if (theEventNr==4 & M1_Size>0) { //d M1
                
                                NR_Deaths++;
                M1_Size--;
                M_Size=M1_Size+M2_Size;
                N=W_Size+M1_Size+M2_Size;
                if (theEnv==1){
                    B_rate_A=BIRTH_RATE_W_E1*W_Size; //*(1-N/K);
                    D_rate_A=DEATH_RATE_W*W_Size;
                    B_rate_M1=BIRTH_RATE_M1_E1*M1_Size; //*(1-N/K);
                    B_rate_M2=BIRTH_RATE_M2_E1*M2_Size; //*(1-N/K);
                    D_rate_M1=DEATH_RATE_M*M1_Size;
                    D_rate_M2=DEATH_RATE_M*M2_Size;
                    SUM_Rates=B_rate_A+D_rate_A+B_rate_M1+B_rate_M2+D_rate_M1+D_rate_M2;
                }
                
                
                
            } else if (theEventNr==5 & M2_Size>0){ //d M2
                
                
                NR_Deaths++;
                M2_Size--;
                M_Size=M1_Size+M2_Size;
                N=W_Size+M1_Size+M2_Size;
                if (theEnv==1){
                    B_rate_A=BIRTH_RATE_W_E1*W_Size; //*(1-N/K);
                    D_rate_A=DEATH_RATE_W*W_Size;
                    B_rate_M1=BIRTH_RATE_M1_E1*M1_Size; //*(1-N/K);
                    B_rate_M2=BIRTH_RATE_M2_E1*M2_Size; //*(1-N/K);
                    D_rate_M1=DEATH_RATE_M*M1_Size;
                    D_rate_M2=DEATH_RATE_M*M2_Size;
                    SUM_Rates=B_rate_A+D_rate_A+B_rate_M1+B_rate_M2+D_rate_M1+D_rate_M2;
                }
            } else { printf("No event works; pop size related;\n"); break; }

        } // end the while.
        // printf("%.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n", B_rate_A, D_rate_A, B_rate_M1, B_rate_M2, D_rate_M1, D_rate_M2);
        timetoextinctall[run]=event_time_total;
        timetoextinctgenerations[run]=NR_Births/(double) K;
        timetoextinctgenerationsdeaths[run]=NR_Deaths/(double) K;
        //record what happened at this simulation
        if (havebreaked==1){
            threshhold++;
        } else if (N==0){
            extinctions++;
        } else {
            printf("wtf?");
        }
    } // End runs of the simulation.
    
    //Write the results to file. -- NOW
    double sumtime=0.0f;
    double meantime=0.0f;
    double sumtimegenerations=0.0f;
    double meantimegenerations=0.0f;
    double sumtimegenerationsdeaths=0.0f;
    double meantimegenerationsdeaths=0.0f;
    
    for (i = 0; i < NRRUNS; i++)
    {
        sumtime = sumtime + timetoextinctall[i];
        sumtimegenerations=sumtimegenerations+timetoextinctgenerations[i];
        sumtimegenerationsdeaths=sumtimegenerationsdeaths+timetoextinctgenerationsdeaths[i];
    }
    meantime=sumtime/NRRUNS;
    meantimegenerations=sumtimegenerations/NRRUNS;
    meantimegenerationsdeaths=sumtimegenerationsdeaths/NRRUNS;
    /*  Compute  variance  and standard deviation  */
    
    double sum1=0.0f;
    double sum2=0.0f;
    double sum3=0.0f;
    for (i = 0; i < NRRUNS; i++)
    {
        sum1 = sum1 + pow((timetoextinctall[i] -meantime), 2);
        sum2 = sum2 + pow((timetoextinctgenerations[i] -meantimegenerations), 2);
        sum3 = sum3 + pow((timetoextinctgenerationsdeaths[i] -meantimegenerationsdeaths), 2);
    }
    double variance = sum1 / (double) NRRUNS;
    double std_deviation = sqrt(variance);
    double variancegenerations = sum2 / (double) NRRUNS;
    double variancegenerationsdeaths = sum3 / (double) NRRUNS;
    double std_deviationgenerations = sqrt(variancegenerations);
    
    double pextinct=0.0f;
    pextinct=extinctions / NRRUNS;
    printf( " ext %d  \n", extinctions);
    printf(" no ext %d  \n", threshhold);
    printf(" no ext freq %f  \n", (float)threshhold/NRRUNS);
    
    char filename[500];
    sprintf(filename, "log/fig2AB_dynamic/outfile_K_%d_memory_%f_BWE1_%f_BM1E1_%f_BM2E1_%f_mutationAa_%f_N1000.txt", K, p, BIRTH_RATE_W_E1, BIRTH_RATE_M1_E1, BIRTH_RATE_M2_E1, mutationAa);
    outfile = fopen(filename,"w");
    
    if (outfile == NULL) {
        fprintf(stderr, "Can't open output file!\n");
        exit(1);
    } else{
        
        fprintf(outfile, "Extinctions %d  \nNoResolution %d \nMeanTime  %f \nVarianceTime  %f \nMeanTimeGenerations %f \nVarianceTimeGenerations %f \nMeanTimeGenerationsDeaths %f \nVarianceTimeGenerationsDeaths %f \n", extinctions, threshhold, meantime, variance, meantimegenerations, variancegenerations,  meantimegenerationsdeaths, variancegenerationsdeaths);
    }
    
    
    // Print simulation total time.
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Time %f\n", time_spent);
    
    gsl_rng_free(gBaseRand); //free rng
    return 0;
    
}
