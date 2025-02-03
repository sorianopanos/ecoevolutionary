#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#define normRAN 2.3283063671E-10F


void inicializarandom();
void copiavector(int *, int **, int);
void monteCarlo_step(int **estados, double **variant_infectivity, double **variant_immunity,int *ninfec, int *nrecup, double mu,double alpha0, double b,int interacciones,double Dinfectivity,double Dimmunity, double mutation_probability, double cross_immunity, int humanostotal, char *);
//void calculagradototal();
void inicializaestados(int **,int *,double rho, int);
unsigned int rueda[256];
unsigned char ind_ran,ig1,ig2,ig3;
void iniciainfectadosrecuperados();
double computeavg(double*,int);
double computevariance(double*,int);
double computeavg_lambda(double*,int*,int size);
double computevariance_lambda(double*,int*,int size);

void iniciahumanos();
void modifyvariant(double, double, int*, double**, double **,int, double,char *);
void inicializavariantes(double **,double**,double,int);
void inicializalambda(double **, double, int);

int main(int kargs,char**argv)
{
    FILE *f,*g;
    int realization,nrealizations,i,j,time,interacciones,humanostotal,Niter,ninfec,nrecup;
    int *estados;
    double x,rectot,lambda0,I0,Dinfectivity,Dimmunity,mu,R0,b,alpha0;
    double *lambda,*variant_immunity,*recovery_probability,variant_avg,lambda_variance,variant_variance,mutation_probability,cross_immunity;
    char aux[600], *flag_model,*flag_execution,*path,output_file[200],output_file2[200];

    inicializarandom();
    //calculagradototal();
    // Arguments got from run files 
    mu=atof(argv[1]); // Infectious period
    Dinfectivity=atof(argv[2]); // Speed evolution of the Infectivity axis
    Dimmunity=atof(argv[3]); // Speed evolution of the Immunity axis
    interacciones=atoi(argv[4]); // Number of contacts (MF approach)
    R0=atof(argv[5]); // Baseline basic reproduction number
    nrealizations=atoi(argv[6]); // Number of realization to get averages
    alpha0=atof(argv[7]); // Baseline virulence rate
    mutation_probability=atof(argv[8]);
    cross_immunity=atof(argv[9]);
    humanostotal=atoi(argv[10]);
    Niter=atoi(argv[11]);
    flag_model=argv[12]; // deterministic/stochastic/tradeoff
    flag_execution=argv[13]; //endemic/trajectories
    path=argv[14];

    I0=0.001;
    estados= malloc(humanostotal*sizeof(int));
    if(strcmp(flag_execution,"endemic")==0)
    {
        strcpy(output_file,path);
        sprintf(aux,"MCMF_%s_%s_mu%.2lf_Dinfectivity%.5lf_Dimmunity%.5lf_k%d_R0%.2lf_alpha0%.3lf_mutationprob%.4lf_crossimmunity%.2lf_N%d.txt",flag_model,flag_execution,mu,Dinfectivity,Dimmunity,interacciones,R0,alpha0,mutation_probability,cross_immunity,humanostotal);
        strcat(output_file,aux);
        f=fopen(output_file,"wt");
    }
    
    if(strcmp(flag_execution,"endemicmorepoints")==0)
    {
        strcpy(output_file,path);
        sprintf(aux,"MCMF_%s_%s_mu%.2lf_Dinfectivity%.5lf_Dimmunity%.5lf_k%d_R0%.2lf_alpha0%.3lf_mutationprob%.4lf_crossimmunity%.2lfN%d.txt",flag_model,flag_execution,mu,Dinfectivity,Dimmunity,interacciones,R0,alpha0,mutation_probability,cross_immunity,humanostotal);
        strcat(output_file,aux);
        f=fopen(output_file,"wt");
    }
    if(strcmp(flag_execution,"trajectories")==0)
    {
       strcpy(output_file,path);
       sprintf(aux,"MCMF_%s_%s_infected_mu%.3lf_Dinfectivity%.5lf_Dimmunity%.5lf_k%d_R0%.2lf_alpha0%.4lf_mutationprob%.4lf_crossimmunity%.2lfN%d.txt",flag_model,flag_execution,mu,Dinfectivity,Dimmunity,interacciones,R0,alpha0,mutation_probability,cross_immunity,humanostotal);
       strcat(output_file,aux);
       f=fopen(output_file,"wt"); 
       fclose(f);    

       strcpy(output_file2,path);
       sprintf(aux,"MCMF_%s_%s_variants_mu%.3lf_Dinfectivity%.5lf_Dimmunity%.5lf_k%d_R0%.2lf_alpha0%.4lf_mutationprob%.4lf_crossimmunity%.2lfN%d.txt",flag_model,flag_execution,mu,Dinfectivity,Dimmunity,interacciones,R0,alpha0,mutation_probability,cross_immunity,humanostotal);
       strcat(output_file2,aux);
       g=fopen(output_file2,"wt");  
       fclose(g);
    }

    if(strcmp(flag_model,"tradeoff")==0 || strcmp(flag_model,"tradeoff2")==0)
    {
        b=R0*(alpha0+mu)/(interacciones*sqrt(alpha0));
        lambda0=b*sqrt(alpha0);
    }
    else
    {
        b=0;
        lambda0=R0*(mu+alpha0)/interacciones;
    }

    for (realization=0;realization<nrealizations;realization++)
    {   
        // Initialize the state
        inicializaestados(&estados,&ninfec,I0,humanostotal);
        inicializavariantes(&lambda,&variant_immunity,lambda0,humanostotal);
        nrecup=0;
        for(time=0;time<Niter;time++)
        {
            if(ninfec!=0)
                monteCarlo_step(&estados,&lambda,&variant_immunity,&ninfec,&nrecup,mu,alpha0,b,interacciones,Dinfectivity,Dimmunity,mutation_probability,cross_immunity,humanostotal,flag_model);
            else
                break;
            
            if(time%5==0)
            {
                if(strcmp(flag_execution,"trajectories")==0)
                {   
                    f=fopen(output_file,"at");
                    g=fopen(output_file2,"at");  
                    fprintf(f,"%d %d %.4lf %.4lf %.4lf\n",realization,time,1.0-(ninfec+nrecup)*1.0/humanostotal,ninfec*1.0/humanostotal,nrecup*1.0/humanostotal);
                    fprintf(g,"%d %d %.4lf %.4lf %.4lf %.4lf\n",realization,time,computeavg_lambda(lambda,estados,humanostotal),computevariance_lambda(lambda,estados,humanostotal),computeavg(variant_immunity,humanostotal),computevariance(variant_immunity,humanostotal));
                    fclose(f);
                    fclose(g);
                }
            }    

            if(time%20==0)
            {
                if(strcmp(flag_execution,"endemicmorepoints")==0)
                {
                    f=fopen(output_file,"at");
                    fprintf(f,"%d %d %.4lf %.4lf %.4lf\n",realization,time,1.0-(ninfec+nrecup)*1.0/humanostotal,ninfec*1.0/humanostotal,nrecup*1.0/humanostotal);
                    fclose(f);
                }

            }
        }
        if(strcmp(flag_execution,"endemic")==0)
        {
            f=fopen(output_file,"at");
            fprintf(f,"%d %d %.4lf %.4lf %.4lf\n",realization,Niter,1.0-(ninfec+nrecup)*1.0/humanostotal,ninfec*1.0/humanostotal,nrecup*1.0/humanostotal);
            fclose(f);
        }
        free(lambda);
        free(variant_immunity);
    }
    return 0;
}

double PariRapu() //devuelve un numero aleatorio entre [0,1) mediante el generador de Parisi-Rapuano
{
    double random;

    ig1=ind_ran-24;
    ig2=ind_ran-55;
    ig3=ind_ran-61;
    rueda[ind_ran]=rueda[ig1]+rueda[ig2];
    random=(rueda[ind_ran]^rueda[ig3]);
    ind_ran++;

    random=random*normRAN;
    return random;
}
void inicializarandom() //inicializa las variables utilizadas en PariRapu()
{
    int i;
    srand(time(NULL));
    ind_ran=ig1=ig2=ig3=0;
    for(i=0;i<256;i++)
        rueda[i]=(rand()<<16)+rand();
}

void copiavector(int* vecoriginal, int **vectorcopia, int humanostotal)
{
    int i,k;
    for(i=0;i<humanostotal;i++)
    {
        (*vectorcopia)[i]=vecoriginal[i];
    }
}

void monteCarlo_step(int **estados,double **lambda,double **variant_immunity,int *ninfec,int *nrecup,double mu,double alpha0, double b,int interacciones,double Dinfectivity,double Dimmunity,double mutation_probability,double cross_immunity,int humanostotal, char *flag_model)
{
    int i,k,j,z,id;
    int controlinfec,numeroinfec,numerorecuperados;
    int *estadosauxh;
    int Niter,finish,control,cont_infect;
    int neighbors[interacciones],idinfec_neighbors[interacciones],ninfec_neighbors;
    double lambda_tot_s,lambda_tot_r,choose_infected_neighbor,auxcum,R0,recovery_probability,prob1;

    estadosauxh=malloc(humanostotal*sizeof(int));
    copiavector(*estados,&estadosauxh,humanostotal);

    for(i=0;i<humanostotal;i++)
    {
        if(strcmp(flag_model,"tradeoff")==0 || strcmp(flag_model,"tradeoff2")==0)
            recovery_probability=mu+((*lambda)[i]*(*lambda)[i])/(b*b);
        else
            recovery_probability=mu+alpha0;

        if((*estados)[i]==1 && PariRapu()<(1-exp(-recovery_probability)))
        {
            if(strcmp(flag_model,"tradeoff2")==0)
            {
                prob1=mu/recovery_probability;
                if(PariRapu()<prob1)
                {
                    estadosauxh[i]=2;
                    (*ninfec)--;
                    (*nrecup)++;
                }
                else
                {
                    estadosauxh[i]=0;
                    (*ninfec)--;
                }
            }
            
            else
            {
                estadosauxh[i]=2;
                (*ninfec)--;
                (*nrecup)++;
            }
        }
        if((*estados)[i]!=1)
        {
            // Get neighbors
            ninfec_neighbors=0;
            lambda_tot_r=0;
            lambda_tot_s=0;
            for(k=0;k<interacciones;k++)
            {
                do
                {
                    id=(int)(PariRapu()*humanostotal);
                } while (id==i);

                if((*estados)[id]==1)
                {
                        idinfec_neighbors[ninfec_neighbors]=id;
                        ninfec_neighbors++;
                        lambda_tot_s+=(*lambda)[id];
                        if((*variant_immunity)[id]>(*variant_immunity)[i])
                            lambda_tot_r+=(*lambda)[id]*(1-exp(-fabs((*variant_immunity)[i]-(*variant_immunity)[id])/cross_immunity));
                }
            }
            // Compute probability of infection for susceptible   
            if((*estados)[i]==0)
            { 
                if(PariRapu()<1-exp(-lambda_tot_s))
                {
                    choose_infected_neighbor=PariRapu()*lambda_tot_s;
                    auxcum=0;
                    for(z=0;z<ninfec_neighbors;z++)
                    {
                            auxcum+=(*lambda)[idinfec_neighbors[z]];
                            if(auxcum>choose_infected_neighbor)
                            {
                                estadosauxh[i]=1;
                                (*variant_immunity)[i]=(*variant_immunity)[idinfec_neighbors[z]];
                                (*lambda)[i]=(*lambda)[idinfec_neighbors[z]];
                                (*ninfec)++;
                                break;
                            }
                    }
                }
            }
            else
            {
                if(PariRapu()<1-exp(-lambda_tot_r))
                {
                    choose_infected_neighbor=PariRapu()*lambda_tot_r;
                    auxcum=0;
                    for(z=0;z<ninfec_neighbors;z++)
                    {
                            if((*variant_immunity)[idinfec_neighbors[z]]>(*variant_immunity)[i])
                                auxcum+=(*lambda)[idinfec_neighbors[z]]*(1-exp(-fabs((*variant_immunity)[i]-(*variant_immunity)[idinfec_neighbors[z]])/cross_immunity));
                            if(auxcum>choose_infected_neighbor)
                            {
                                estadosauxh[i]=1;
                                (*variant_immunity)[i]=(*variant_immunity)[idinfec_neighbors[z]];
                                (*lambda)[i]=(*lambda)[idinfec_neighbors[z]];
                                (*ninfec)++;
                                (*nrecup)--;
                                break;
                            }
                    }
                }
            }
        }
    }
        copiavector(estadosauxh,estados,humanostotal);
        modifyvariant(Dinfectivity,Dimmunity,*estados,lambda,variant_immunity,humanostotal,mutation_probability,flag_model);
        free(estadosauxh);       
}
void inicializaestados(int ** estados,int *ninfec,double rho,int humanostotal)
{

    int j;
    *ninfec=0;
    for(j=0;j<humanostotal;j++){
        if(PariRapu()<rho)
        {
            (*estados)[j]=1;
            (*ninfec)++;           
        
        }
        else
            (*estados)[j]=0;

    }
}


void modifyvariant(double Dinfectivity, double Dimmunity, int *estados, double **variant_infectivity, double **variant_immunity,int humanostotal,double mutation_probability, char *flag_model)
{
  int i;
  double r1,r2,r3,r4;
  for (i=0;i<humanostotal;i++)
  {
    if(estados[i]==1)
    {
        if (PariRapu()<mutation_probability)
        {
            if(strcmp(flag_model,"deterministic")!=0)
            {
                r1=r2=0;
                do{
                r1=PariRapu();
                }
                while(r1==0);

                do{
                r2=PariRapu();
                }
                while(r2==0);

                do{
                r3=PariRapu();
                }
                while(r3==0);

                do{
                r4=PariRapu();
                }
                while(r4==0);

                (*variant_immunity)[i]=(*variant_immunity)[i]+Dimmunity*sqrt(-2*log(r1))*cos(2*M_PI*r2);
                (*variant_infectivity)[i]=(*variant_infectivity)[i]+Dinfectivity*sqrt(-2*log(r3))*cos(2*M_PI*r4);
                if((*variant_infectivity)[i]<0)
                    (*variant_infectivity)[i]=0;
            }
            else
            {
                (*variant_immunity)[i]=(*variant_immunity)[i]+Dimmunity;
                (*variant_infectivity)[i]=(*variant_infectivity)[i]+Dinfectivity;
            }
        }
    }
  }
}

void inicializavariantes(double **lambda, double **variant_immunity,double lambda0,int humanostotal)
{
    int i;

    (*lambda)=malloc(humanostotal*sizeof(double));
    (*variant_immunity)=malloc(humanostotal*sizeof(double));
    for(i=0;i<humanostotal;i++)
    {
        (*lambda)[i]=lambda0;
        (*variant_immunity)[i]=0;
    }
}


double computeavg(double *vector,int size)
{
    int i;
    double mean;
    mean=0;
    for(i=0;i<size;i++)
     mean+=vector[i];

    mean/=size;
    return mean;
}

double computevariance(double *vector,int size)
{
    int i;
    double mean,mean_sq,variance;
    mean=computeavg(vector,size);
    mean_sq=0;
    for(i=0;i<size;i++)
     mean_sq+=vector[i]*vector[i];

    variance=mean_sq/size-mean*mean;

    return variance;
}


double computeavg_lambda(double *vector,int *estados,int size)
{
    int i,size_mean;
    double mean;
    mean=0;
    size_mean=0;
    for(i=0;i<size;i++)
        if(estados[i]==1){   
            mean+=vector[i];
            size_mean++;
            }
    if(size_mean!=0)
        mean/=size_mean;

    return mean;
}

double computevariance_lambda(double *vector,int *estados,int size)
{
    int i,size_mean;
    double mean,mean_sq,variance;
    mean=computeavg_lambda(vector,estados,size);
    mean_sq=0;
    size_mean=0;
    for(i=0;i<size;i++)
        if(estados[i]==1)
        {
            mean_sq+=vector[i]*vector[i];
            size_mean++;
        }
    variance=mean_sq/size_mean-mean*mean;

    return variance;
}


